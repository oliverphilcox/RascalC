"""
Function to post-process the single-field 3PCF Legendre binned integrals computed by the C++ code, obtaining the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.
"""

import numpy as np
import numpy.typing as npt
from scipy.optimize import fmin
import os
from ..utils import format_skip_r_bins
from ..post_process.utils import check_eigval_convergence, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_3pcf_legendre
from .utils import cov_filter_3pcf_legendre, load_matrices, add_cov_terms
from typing import Callable, Iterable


def Psi(alpha: float, c3: npt.NDArray[np.float64], c4: npt.NDArray[np.float64], c5: npt.NDArray[np.float64], c6: npt.NDArray[np.float64], c3s: npt.NDArray[np.float64], c4s: npt.NDArray[np.float64], c5s: npt.NDArray[np.float64], c6s: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
    c_tot = add_cov_terms(c3, c4, c5, c6, alpha)
    partial_covs = add_cov_terms(c3s, c4s, c5s, c6s, alpha)
    _, Psi = compute_D_precision_matrix(partial_covs, c_tot)
    return Psi


def neg_log_L1(alpha: float, target_cov: npt.NDArray[np.float64], c3: npt.NDArray[np.float64], c4: npt.NDArray[np.float64], c5: npt.NDArray[np.float64], c6: npt.NDArray[np.float64], c3s: npt.NDArray[np.float64], c4s: npt.NDArray[np.float64], c5s: npt.NDArray[np.float64], c6s: npt.NDArray[np.float64]) -> float:
    """Return negative log L1 likelihood between 3PCF theory and target (data jackknife or mock sample) covariance matrices.
    log L1 is the Kullback-Leibler divergence with constant terms (including log(det(target_cov))) removed.
    As a result, the `target_cov` can be a singular matrix.
    This function does not allow negative shot-noise rescaling `alpha` by returning infinity."""
    if alpha < 0: return np.inf # negative shot-noise rescaling causes problems and does not make sense
    Psi_alpha = Psi(alpha, c3, c4, c5, c6, c3s, c4s, c5s, c6s)
    logdet = np.linalg.slogdet(Psi_alpha)
    if logdet[0] < 0:
        # Remove any dodgy inversions
        return np.inf        
    return np.trace(np.matmul(Psi_alpha, target_cov)) - logdet[1]


def fit_shot_noise_rescaling(target_cov: npt.NDArray[np.float64], c3: npt.NDArray[np.float64], c4: npt.NDArray[np.float64], c5: npt.NDArray[np.float64], c6: npt.NDArray[np.float64], c3s: npt.NDArray[np.float64], c4s: npt.NDArray[np.float64], c5s: npt.NDArray[np.float64], c6s: npt.NDArray[np.float64]) -> float:
    """Fit the 3PCF covariance matrix model to `target_cov` to find the optimal shot-noise rescaling.
    `target_cov` can be a singular matrix."""
    alpha_best = fmin(neg_log_L1, 1., args = (target_cov, c3, c4, c5, c6, c3s, c4s, c5s, c6s))
    return alpha_best[0]


def post_process_3pcf_mocks(mock_cov_file: str, file_root: str, n: int, max_l: int, outdir: str | None = None, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, exclude_samebins: bool = True, exclude_odd_l: bool = False, check_finished: bool = True, max_l_mock: int | None = None, print_function: Callable[[str], None] = print, dry_run: bool = False) -> dict[str]:
    r"""
    3PCF post-processing for Legendre (accumulated) mode, obtaining the shot-noise rescaling parameter, alpha, from a mock-derived covariance matrix.

    Now it should be safe to run this post-processing while the main RascalC computation is still running, as long as you do not put multiple runs into one output directory.
    This is achieved by a default heuristic check for the normal finishing of the main RascalC computation.
    With this, inspecting the output of an aborted or timed-out run is harder, by default it will be considered unfinished. But if you are sure that the main computation is not running, you can disable the check via ``check_finished=False``. Doing this once should be sufficient, repeated post-processing attempts should no longer detect the run as unfinished.

    Parameters
    ----------
    mock_cov_file : string
        Path to the text file containing the mock sample covariance matrix (from ENCORE 3PCF measurements), which will be used to determine the shot-noise rescaling parameter, alpha. The covariance matrix should have the same radial bins as the theoretical matrices computed by RascalC, have the repeated radial bin pairs removed and same-bin pairs removed according to the ``exclude_samebins`` option (typically, yes). It may have a different (higher) maximum ell than the computed ``RascalC`` covariance matrices (in such a case, specify it with max_l_mock). The ``skip_r_bins`` and ``skip_l`` filters will be applied after the mock covariance matrix is loaded. The covariance matrix should be in a text format that can be loaded by ``numpy.loadtxt``.
    
    file_root : string
        Path to the RascalC (:func:`RascalC.run_cov_3pcf` or command-line) output directory.
    
    n : integer
        The number of radial bins used in the RascalC run (before applying ``skip_r_bins`` if it is provided).
    
    max_l : integer
        The maximum ell (Legendre moment index) used in the RascalC run (before applying ``skip_l`` if it is provided).

    outdir : string or None
        (Optional) path to the directory in which the post-processing results should be saved. If None (default), is set to ``file_root``. Empty string means the current working directory.
        We advise to use different output directories for different post-processing options.

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
        (Optional) If False, the covariance will include the pairs of the same radial bins. This probably will not work with the mock sample covariance file produced from ENCORE.
        The default behavior (for the True value) is to exclude them for compatibility with ENCORE.
        In either case, the post-processed covariances only include each pair of different radial bins in one ordering, ``bin1 < bin2``; the raw covariances also include ``bin1 > bin2`` pairs.
    
    exclude_odd_l : boolean
        (Optional) If True, the covariance will exclude the odd multipoles; note that then they will also not count in ``skip_l``. By default (False value), odd multipoles are kept and counted in ``skip_l``.
    
    max_l_mock : integer or None
        (Optional) The maximum ell (Legendre moment index) used in the mock covariance matrix. If None (default), it is assumed to be the same as `max_l` for theoretical covariance matrices. This option is provided in case the mock covariance matrix has a different (higher) maximum ell than the theoretical covariance matrices computed by RascalC.
    
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

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_3PCF_Mocks_n%d_l%d.npz' % (n, max_l))
    name_dict = dict(path=output_name, filename=os.path.basename(output_name))
    if dry_run: return name_dict

    # Load the mock covariance matrix and convert it to the same convention as RascalC theoretical covariance matrices
    mock_cov = np.loadtxt(mock_cov_file) # load external mock covariance matrix from text file, would be based on ENCORE
    # prepare r bin filtering
    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    r_bin_indices = np.arange(n)
    r_bin_index1 = np.repeat(r_bin_indices, n-exclude_samebins-r_bin_indices)
    r_bin_index2 = np.concatenate([r_bin_indices[i+exclude_samebins:] for i in r_bin_indices])
    n_r_pairs_orig = len(r_bin_index1)
    # r_bin_index1 and r_bin_index2 cover all the bin pairs under the condition r_bin_index1 < r_bin_index2, the order follows the ENCORE format
    r_filter = (r_bin_index1 >= skip_r_bins_start) & (r_bin_index1 < n - skip_r_bins_end) & (r_bin_index2 >= skip_r_bins_start) & (r_bin_index2 < n - skip_r_bins_end) # filter for the bin pairs to keep based on the skip_r_bins option
    n_r_pairs = r_filter.sum()
    # prepare ell indexing and scaling factor accounting for the different basis
    ells = np.arange(0, max_l + 1, 1+exclude_odd_l)
    if skip_l > 0: ells = ells[:-skip_l] # without the condition, wouldn't work right for skip_l=0
    ell_factor = ((-1)**ells * np.sqrt(2 * ells + 1) / (4 * np.pi)) # the ell-dependent factor between the ENCORE 3-point basis functions and Legendre polynomials given by Equation (17) in https://arxiv.org/pdf/2105.08722
    # scale and transpose the covariance
    if max_l_mock is None: max_l_mock = max_l
    mock_cov = mock_cov.reshape(max_l_mock + 1, n_r_pairs_orig, max_l_mock + 1, n_r_pairs_orig) # reshape the covariance from 2D to 4D, the ENCORE ordering is [l, r_bin_pair] for both rows and columns
    mock_cov = mock_cov[ells][:, :, ells] # apply the multipole selection
    mock_cov = mock_cov.transpose(1, 0, 3, 2) # change ordering to [r_bin_pair, l] for both rows and columns as in RascalC
    mock_cov *= ell_factor[:, None, None] * ell_factor[None, None, :] # apply the factor, NumPy broadcasting matches trailing dimensions. need to double-check if it is not division; there might also be a factor of 2 or something similar
    mock_cov = mock_cov[r_filter][:, :, r_filter] # apply the r bin pair filter
    mock_cov = mock_cov.reshape(n_r_pairs * len(ells), n_r_pairs * len(ells)) # reshape back to 2D, should be ready

    cov_filter = cov_filter_3pcf_legendre(n, max_l, skip_r_bins, skip_l, exclude_samebins, exclude_odd_l)
    
    input_file = load_raw_covariances_3pcf_legendre(file_root, n, max_l, n_samples, check_finished=check_finished, print_function=print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c3, c4, c5, c6 = load_matrices(input_file, n, max_l, cov_filter, full=True)

    # Check matrix convergence by analogy with 2PCF, may be less helpful
    eigval_ok = check_eigval_convergence(c3, c6, Npcf=3, print_function=print_function)

    # Load in partial theoretical matrices
    c3s, c4s, c5s, c6s = load_matrices(input_file, n, max_l, cov_filter, full=False)

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    alpha_best = fit_shot_noise_rescaling(mock_cov, c3, c4, c5, c6, c3s, c4s, c5s, c6s)
    print_function("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)

    # Check matrix convergence for the optimal alpha: if it is <1, the eigenvalue criterion should be strengthened
    if eigval_ok and alpha_best < 1: check_eigval_convergence(c3, c6, alpha_best, Npcf=3, print_function=print_function)

    # Compute full covariance matrix
    full_cov = add_cov_terms(c3, c4, c5, c6, alpha_best)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    partial_cov = add_cov_terms(c3s, c4s, c5s, c6s, alpha_best)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)  

    output_dict = dict(full_theory_covariance=full_cov, shot_noise_rescaling=alpha_best,
                       full_theory_precision=full_prec, N_eff=N_eff_D,
                       full_theory_D_matrix=full_D_est, individual_theory_covariances=partial_cov,
                       mock_covariance=mock_cov)
    
    np.savez_compressed(output_name, **output_dict)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict | name_dict