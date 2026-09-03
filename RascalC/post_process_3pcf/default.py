"""
Function to post-process the single-field 3PCF Legendre binned integrals computed by the C++ code with a given shot-noise rescaling parameter value, alpha.
We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.
"""

import numpy as np
import os
from ..post_process.utils import check_eigval_convergence, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_3pcf_legendre
from .utils import cov_filter_3pcf_legendre, load_matrices, add_cov_terms
from typing import Callable, Iterable


def post_process_3pcf(file_root: str, n: int, max_l: int, outdir: str | None = None, alpha: float = 1, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, exclude_samebins: bool = True, exclude_odd_l: bool = False, check_finished: bool = True, print_function: Callable[[str], None] = print, dry_run: bool = False) -> dict[str]:
    r"""
    3PCF post-processing for Legendre (accumulated) mode for a given shot-noise rescaling parameter value, alpha.

    Now it should be safe to run this post-processing while the main RascalC computation is still running, as long as you do not put multiple runs into one output directory.
    This is achieved by a default heuristic check for the normal finishing of the main RascalC computation.
    With this, inspecting the output of an aborted or timed-out run is harder, by default it will be considered unfinished. But if you are sure that the main computation is not running, you can disable the check via ``check_finished=False``. Doing this once should be sufficient, repeated post-processing attempts should no longer detect the run as unfinished.

    Parameters
    ----------
    file_root : string
        Path to the RascalC (:func:`RascalC.run_cov_3pcf` or command-line) output directory.
    
    n : integer
        The number of radial bins used in the RascalC run (before applying ``skip_r_bins`` if it is provided).
    
    max_l : integer
        The maximum ell (Legendre moment index) used in the RascalC run (before applying ``skip_l`` if it is provided).

    outdir : string or None
        (Optional) path to the directory in which the post-processing results should be saved. If None (default), is set to ``file_root``. Empty string means the current working directory.
        We advise to use different output directories for different post-processing options.

    alpha : float
        Fixed shot-noise rescaling value to use. In principle optional, but the default value of 1 may not be particularly good.

    skip_r_bins : integer or tuple of two integers
        (Optional) removal of some radial bins.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.

    skip_l : integer
        (Optional) number of higher multipoles to skip (from the end, counting all multipoles by default and only even multipoles if exlude_odd_l is True).

    n_samples : None, integer, array/list/tuple/etc of integers or boolean values
        (Optional) selection of RascalC subsamples (independent realizations of Monte-Carlo integrals).
        
            - If None, use all (default).
            - If an integer, use the given number of samples from the beginning.
            - If an array/list/tuple/etc of integers, it will be used as a NumPy index array.
            - If an array/list/tuple/etc of boolean, it will be used as a NumPy boolean array mask.
    
    exclude_samebins : boolean
        (Optional) If False, the covariance will include the pairs of the same radial bins.
        The default behavior (for the True value) is to exclude them for compatibility with ENCORE.
        In either case, the post-processed covariances only include each pair of different radial bins in one ordering, ``bin1 < bin2``; the raw covariances also include ``bin1 > bin2`` pairs.
    
    exclude_odd_l : boolean
        (Optional) If True, the covariance will exclude the odd multipoles; note that then they will also not count in ``skip_l``. By default (False value), odd multipoles are kept and counted in ``skip_l``.
    
    print_function : Callable
        (Optional) custom function to use for printing. Default is ``print``.
    
    dry_run: boolean
        (Optional) If True, this will not run actual post-processing, only determine the filename and path (see below).

    Returns
    -------
    post_processing_results : dict[str]
        Post-processing results as a dictionary with string keys and Numpy array values (mostly). All this information is also saved in a ``Rescaled_Covariance_Matrices*.npz`` file in the ``out_dir`` (in ``file_root`` if the former is not provided).
        Selected common keys are: ``"full_theory_covariance"`` for the final covariance matrix and ``"shot_noise_rescaling"`` for the shot-noise rescaling value(s).
        For convenience, in the output dictionary only, ``"filename"`` contains the name of the file where the results were saved (which can be inconvenient to predict), and ``"path"`` contains its path (also obtainable by :func:`os.path.join`-ing ``out_dir`` with the filename)
    """
    # Set default output directory if not set
    if outdir is None: outdir = file_root

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_3PCF_n%d_l%d.npz' % (n, max_l))
    name_dict = dict(path=output_name, filename=os.path.basename(output_name))
    if dry_run: return name_dict

    cov_filter = cov_filter_3pcf_legendre(n, max_l, skip_r_bins, skip_l, exclude_samebins, exclude_odd_l)
    
    input_file = load_raw_covariances_3pcf_legendre(file_root, n, max_l, n_samples, check_finished=check_finished, print_function=print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c3, c4, c5, c6 = load_matrices(input_file, n, max_l, cov_filter, full=True)

    # Check matrix convergence by analogy with 2PCF, may be less helpful
    check_eigval_convergence(c3, c6, alpha, Npcf=3, print_function=print_function)

    # Compute full covariance matrix
    full_cov = add_cov_terms(c3, c4, c5, c6, alpha)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial theoretical matrices
    c3s, c4s, c5s, c6s = load_matrices(input_file, n, max_l, cov_filter, full=False)
    partial_cov = add_cov_terms(c3s, c4s, c5s, c6s, alpha)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)  

    output_dict = dict(full_theory_covariance=full_cov, shot_noise_rescaling=alpha,
                       full_theory_precision=full_prec, N_eff=N_eff_D,
                       full_theory_D_matrix=full_D_est, individual_theory_covariances=partial_cov)
    
    np.savez_compressed(output_name, **output_dict)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict | name_dict