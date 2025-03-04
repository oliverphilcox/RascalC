from typing import Iterable
import os
from glob import glob
from re import fullmatch
from warnings import warn
from ..interface import indices_corr_all
from ..raw_covariance_matrices import collect_raw_covariance_matrices
from ..convergence_check_extra import convergence_check_extra
from .default import post_process_default
from .default_multi import post_process_default_multi
from .jackknife import post_process_jackknife
from .jackknife_multi import post_process_jackknife_multi
from .legendre import post_process_legendre
from .legendre_multi import post_process_legendre_multi
from .legendre_mix_jackknife import post_process_legendre_mix_jackknife


def post_process_auto(file_root: str, out_dir: str | None = None, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, tracer: int = 1, n_samples: None | int | Iterable[int] | Iterable[bool] = None, shot_noise_rescaling1: float = 1, shot_noise_rescaling2: float = 1, print_function = print, extra_convergence_check: bool = True, jackknife: bool | None = None, legendre: bool | None = None, two_tracers: bool | None = None, n_r_bins: int | None = None, n_mu_bins: int | None = None, n_jack: int | None = None, max_l: int | None = None) -> dict[str]:
    """
    Automatic but highly customizable post-processing interface. Designed to work with the ``RascalC.run_cov`` outputs.

    Parameters
    ----------
    file_root: string
        Path to the RascalC output directory.

    out_dir: string | None
        (Optional) path to the directory in which the post-processing results should be saved. If None (default), is set to ``file_root``. Empty string means the current working directory.

    skip_r_bins: integer or tuple of two integers
        (Optional) removal of some radial bins.
        First number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.

    skip_l: integer
        (Optional) number of higher multipoles to skip (from the end).

    tracer: 1 or 2
        (Optional) if the RascalC output directory contains two-tracer results, this allows to select the second tracer for single-tracer post-processing.

    n_samples: None, integer, array/list/tuple/etc of integers or boolean values
        (Optional) selection of RascalC subsamples (independent realizations of Monte-Carlo integrals).
        If None, use all (default).
        If an integer, use the given number of samples from the beginning.
        If an array/list/tuple/etc of integers, it will be used as a NumPy index array.
        If an array/list/tuple/etc of boolean, it will be used as a NumPy boolean array mask.

    shot_noise_rescaling1: float
        (Optional) shot-noise rescaling value for the first tracer (default 1). Auto-determined with jackknife.

    shot_noise_rescaling2: float
        (Optional) shot-noise rescaling value for the second tracer only in multi-tracer mode (default 1). Auto-determined with jackknife.
    
    print_function: Callable
        (Optional) custom function to use for printing. Default is ``print``.

    extra_convergence_check: bool
        (Optional) whether to perform the extra convergence check. It is done by default.

    jackknife: boolean or None
        (Optional) boolean value sets jackknife mode manually. If None (default), this mode is determined automatically.

    legendre: boolean or None
        (Optional) boolean value sets Legendre mode manually. If None (default), this mode is determined automatically.

    two_tracers: boolean or None
        (Optional) boolean value sets 1- vs 2-tracer mode manually. If None (default), this mode is determined automatically.
    
    n_r_bins, n_mu_bins, n_jack, max_l: integer or None
        (Optional) integer value is used to set manually the number of radial bins, angular bins (not needed in some Legendre modes), jackknife regions (jackknife mode only) and maximum (even) ell (Legendre mode only) respectively. Each parameter which is None (default) is determined automatically.

    Returns
    -------
    post_processing_results : dict[str, np.ndarray[float]]
        Post-processing results as a dictionary with string keys and Numpy array values. All this information is also saved in a ``Rescaled_Covariance_Matrices*.npz`` file in the ``out_dir`` (in ``file_root`` if the former is not provided).
        Selected common keys are: ``"full_theory_covariance"`` for the final covariance matrix and ``"shot_noise_rescaling"`` for the shot-noise rescaling value(s).
    """
    # Set default output directory if not set
    if out_dir is None: out_dir = file_root

    # Simple auto-determination of modes
    if jackknife is None: jackknife = os.path.isdir(os.path.join(file_root, "xi_jack"))
    if two_tracers is None: two_tracers = os.path.isfile(os.path.join(file_root, f"xi/xi_22.dat"))

    ntracers = 2 if two_tracers else 1
    ncorr = ntracers * (ntracers + 1) // 2
    indices_corr = indices_corr_all[:ncorr]

    legendre_orig = len(glob(os.path.join(file_root, f"BinCorrectionFactor*"))) > 0
    legendre_mix = len(glob(os.path.join(file_root, "weights/mu_bin_legendre_factors_*.txt"))) > 0
    if legendre is None: legendre = legendre_orig or legendre_mix

    print_function(f"Legendre: {legendre}")
    print_function(f"Jackknife: {jackknife}")
    print_function(f"Number of tracers: {1 + two_tracers}")

    if legendre_orig and jackknife: warn("Direct accumulation Legendre mode is not compatible with jackknives")

    # Determine number of radial, mu bins and/or jackknives automatically as needed
    binned_pair_names = glob("binned_pair_counts_n*_m*_j*_??.dat" if jackknife else "RR_counts_n*_m*_??.dat", root_dir = os.path.join(file_root, "weights"))
    if len(binned_pair_names) < ncorr: raise ValueError(f"Need {ncorr} pair counts, found {len(binned_pair_names)}")
    rstr = r'binned_pair_counts_n(?P<N_R_BINS>\d+)_m(?P<N_MU_BINS>\d+)_j(?P<N_JACK>\d+)_(?P<CORR_INDEX>\d+).dat' if jackknife else r'RR_counts_n(?P<N_R_BINS>\d+)_m(?P<N_MU_BINS>\d+)_(?P<CORR_INDEX>\d+).dat' # regex
    if (m := fullmatch(rstr, binned_pair_names[0])):
        if n_r_bins is None: n_r_bins = int(m["N_R_BINS"])
        if n_mu_bins is None: n_mu_bins = int(m["N_MU_BINS"])
        if jackknife and n_jack is None: n_jack = int(m["N_JACK"])
    else: warn("The pair count names not matched to the pattern. Not able to autodetermine `n_r_bins`, `n_mu_bins` and `n_jack`.")

    # Determine max_l automatically if needed
    if legendre and max_l is None:
        raw_filenames = glob(f"Raw_Covariance_Matrices_n{n_r_bins}_l*.npz", root_dir = file_root)
        if raw_filenames:
            if len(raw_filenames) > 1: warn("Found multiple `max_l` options.")
            rstr = fr"Raw_Covariance_Matrices_n{n_r_bins}_l(?P<MAX_L>\d+).npz"
            if not (m := fullmatch(rstr, raw_filenames[0])): raise ValueError("Raw covariance matrices filename suddenly not matched")
            max_l = int(m["MAX_L"])
        else:
            prefix = f"n{n_r_bins}_l"
            raw_matrices = collect_raw_covariance_matrices(file_root, print_function = print_function)
            matched_labels = [label[len(prefix):] for label in raw_matrices.keys() if label.startswith(prefix)]
            if not matched_labels: raise ValueError("No Legendre results matched by the number of radial bins.")
            if len(matched_labels) > 1: warn("Found multiple `max_l` options.")
            max_l = int(matched_labels[0])

    if jackknife:
        xi_jack_names = [os.path.join(file_root, f"xi_jack/xi_jack_n{n_r_bins}_m{n_mu_bins}_j{n_jack}_{index}.dat") for index in indices_corr]

    if two_tracers:
        if legendre:
            results = post_process_legendre_multi(file_root, n_r_bins, max_l, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_r_bins, skip_l, n_samples = n_samples, print_function = print_function)
        elif jackknife:
            results = post_process_jackknife_multi(*xi_jack_names, os.path.join(file_root, "weights"), file_root, n_mu_bins, out_dir, skip_r_bins, n_samples = n_samples, print_function = print_function)
        else: # default
            results = post_process_default_multi(file_root, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_r_bins, n_samples = n_samples, print_function = print_function)
    else:
        if legendre:
            if jackknife:
                results = post_process_legendre_mix_jackknife(xi_jack_names[0], os.path.join(file_root, "weights"), file_root, n_mu_bins, max_l, out_dir, skip_r_bins, skip_l, tracer = tracer, n_samples = n_samples, print_function = print_function)
            else:
                results = post_process_legendre(file_root, n_r_bins, max_l, out_dir, shot_noise_rescaling1, skip_r_bins, skip_l, tracer = tracer, n_samples = n_samples, print_function = print_function)
        elif jackknife:
            results = post_process_jackknife(xi_jack_names[0], os.path.join(file_root, "weights"), file_root, n_mu_bins, out_dir, skip_r_bins, tracer = tracer, n_samples = n_samples, print_function = print_function)
        else: # default
            results = post_process_default(file_root, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, skip_r_bins, tracer = tracer, n_samples = n_samples, print_function = print_function)

    if extra_convergence_check:
        print_function("Performing an extra convergence check")
        convergence_check_extra(results, print_function = print_function)

    return results
    