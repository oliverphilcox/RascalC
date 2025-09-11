from typing import Iterable, Callable, Literal
import numpy as np
import os
from glob import glob
from re import fullmatch
from warnings import warn
from ..utils import indices_corr_all
from ..raw_covariance_matrices import collect_raw_covariance_matrices
from ..convergence_check_extra import convergence_check_extra
from .default import post_process_default
from .default_multi import post_process_default_multi
from .jackknife import post_process_jackknife
from .jackknife_multi import post_process_jackknife_multi
from .legendre import post_process_legendre
from .legendre_multi import post_process_legendre_multi
from .legendre_mix_jackknife import post_process_legendre_mix_jackknife
from .legendre_mix_jackknife_multi import post_process_legendre_mix_jackknife_multi

# dependencies for mock post-processing if we want to separate them at some point
import pycorr
from .default_mocks import post_process_default_mocks
from .default_mocks_multi import post_process_default_mocks_multi
from .legendre_mocks import post_process_legendre_mocks
from .legendre_mocks_multi import post_process_legendre_mocks_multi
from ..pycorr_utils.sample_cov_multipoles import sample_cov_multipoles_from_pycorr_to_file
from ..pycorr_utils.sample_cov import sample_cov_from_pycorr_to_file


def post_process_auto(file_root: str,
                      out_dir: str | None = None,
                      skip_s_bins: int | tuple[int, int] = 0, skip_l: int = 0,
                      tracer: Literal[1, 2] = 1,
                      n_samples: None | int | Iterable[int] | Iterable[bool] = None,
                      shot_noise_rescaling1: float = 1, shot_noise_rescaling2: float = 1,
                      xi_11_samples: Iterable[pycorr.twopoint_estimator.BaseTwoPointEstimator] | None = None,
                      xi_12_samples: Iterable[pycorr.twopoint_estimator.BaseTwoPointEstimator] | None = None,
                      xi_22_samples: Iterable[pycorr.twopoint_estimator.BaseTwoPointEstimator] | None = None,
                      xi_sample_cov: np.ndarray[float] | None = None,
                      print_function: Callable[[str], None] = print,
                      extra_convergence_check: bool = True,
                      jackknife: bool | None = None, load_sample_cov: bool | None = None,
                      legendre: bool | None = None, two_tracers: bool | None = None,
                      n_r_bins: int | None = None, n_mu_bins: int | None = None,
                      n_jack: int | None = None,
                      max_l: int | None = None,
                      dry_run: bool = False) -> dict[str]:
    r"""
    Automatic but highly customizable post-processing interface. Designed to work with the :func:`RascalC.run_cov` outputs.

        - By default, this function guesses jackknife or mock pipeline and the covariance binning mode by the output directory contents, but you can also specify some or all of these regimes via optional arguments.
        - Note that ``skip_s_bins`` and ``skip_l`` are not auto-determined, so by default no bins are be skipped even if they were cut in :func:`RascalC.run_cov`.

    Do not run this (or any other post-processing function/script) while the main RascalC computation is running — this may delete the output directory and cause the code to crash.

    Parameters
    ----------
    file_root : string
        Path to the RascalC (:func:`RascalC.run_cov`) output directory.
        This is the only necessary argument. If no others are provided, jackknife or mock pipeline and the covariance binning mode will be determined automatically. The result should match the :func:`RascalC.run_cov` setup (except non-zero ``skip_s_bins`` and ``skip_l``, which are not auto-determined).

    out_dir : string | None
        (Optional) path to the directory in which the post-processing results should be saved. If None (default), is set to ``file_root``. Empty string means the current working directory.
        We advise to use different output directories for different post-processing options.

    skip_s_bins : integer or tuple of two integers
        (Optional) removal of some radial bins.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.

    skip_l : integer
        (Optional) number of higher multipoles to skip (from the end and counting only even multipoles).

    tracer : 1 or 2
        (Optional) if the RascalC output directory contains two-tracer results, ``tracer = 2`` together with ``two_tracers = False`` allows to select the second tracer for single-tracer post-processing.

    n_samples : None, integer, array/list/tuple/etc of integers or boolean values
        (Optional) selection of RascalC subsamples (independent realizations of Monte-Carlo integrals).
        
            - If None, use all (default).
            - If an integer, use the given number of samples from the beginning.
            - If an array/list/tuple/etc of integers, it will be used as a NumPy index array.
            - If an array/list/tuple/etc of boolean, it will be used as a NumPy boolean array mask.

    shot_noise_rescaling1 : float
        (Optional) shot-noise rescaling value for the first tracer (default 1).
        In jackknife or mock mode, the shot-noise rescaling value is auto-determined, so this parameter has no effect.

    shot_noise_rescaling2 : float
        (Optional) shot-noise rescaling value for the second tracer only in multi-tracer mode (default 1).
        In jackknife or mock mode, the shot-noise rescaling values are auto-determined, so this parameter has no effect.
    
    xi_11_samples: None, or list or tuple of ``pycorr.TwoPointEstimator``\s
        (Optional) A set of ``pycorr.TwoPointEstimator``\s (typically from mocks) for the first tracer auto-correlation function to compute the sample covariance to use as reference in shot-noise rescaling optimization.
        Must have the same binning as the covariance (except possibly the angular/mu bins in Legendre mode).
        Providing this option is not compatible with jackknife (enabled explicitly).
        For two-tracer post-processing, providing this requires also passing ``xi_22_samples`` with the same number of samples and binning.
    
    xi_22_samples: None, or list or tuple of ``pycorr.TwoPointEstimator``\s
        (Optional, only for two-tracer post-processing) A set of ``pycorr.TwoPointEstimator``\s (typically from mocks) for the second tracer auto-correlation function to compute the sample covariance to use as reference in shot-noise rescaling optimization.
        For two-tracer post-processing, this is necessary if ``xi_11_samples`` are provided, and the binning and the number of samples must be the same.
    
    xi_12_samples: None, or list or tuple of ``pycorr.TwoPointEstimator``\s
        (Optional, only for two-tracer post-processing) A set of ``pycorr.TwoPointEstimator``\s (typically from mocks) for the cross-correlation function between the first and the second tracers to compute the sample covariance.
        Must have the same binning and number of samples as ``xi_11_samples`` if provided.
        Only auto-covariances of the auto-correlation functions are used in two-tracer mock-based post-processing, so cross-correlation computations may be omitted if not needed for other reasons.
    
    xi_sample_cov: None, or a symmetric positive definite matrix
        (Optional) The pre-computed sample covariance to use as reference in shot-noise rescaling optimization.
        Providing this option is not compatible with ``xi_11_samples`` and jackknife (enabled explicitly).
        Please ensure the right bin ordering if you use this option (because of this, it may be easier to use ``xi_11_samples`` instead).
        In "s_mu" ``mode``, the top-level ordering/grouping is by separation/radial bins and then by angular/µ bins, i.e. the neighboring angular/µ bins (after wrapping!) in one separation/radial bin are next to each other.
        In any of the Legendre ``mode``\s, the top-level ordering/grouping is by multipoles and then by separation/radial bins, i.e. the same multipole moments in neighboring radial bins are next to each other.
        For multi-tracer, the topmost-level ordering must be by the correlation function: 11, 12, 22.
        If you use the ``skip_s_bins`` option (and/or ``skip_l`` with Legendre binning), the requested bins will be removed from the sample covariance. I.e., the computed sample covariance must match the ``RascalC`` binning before those cuts.
    
    print_function : Callable
        (Optional) custom function to use for printing. Default is ``print``.

    extra_convergence_check : bool
        (Optional) whether to perform the extra convergence check. It is done by default.

    jackknife : boolean or None
        (Optional) boolean value sets jackknife mode manually. If None (default), this mode is determined automatically.

    load_sample_cov : boolean or None
        (Optional) boolean value sets whether to load the (mock) sample covariance saved in a default file under ``file_root``. If None (default), this is determined automatically by the existence of the file. Enabling this option contradicts with ``xi_11_samples``, ``xi_sample_cov`` and ``jackknife``.
        Note that auto-determined jackknife disables the automatic mock post-processing (because jackknife pipeline has higher priority). In practice, if your ``file_root`` has both jackknife files and the (mock) sample covariance, you need ``load_sample_cov=True`` to load the latter.

    legendre : boolean or None
        (Optional) boolean value sets Legendre (vs s,mu) mode manually. If None (default), this mode is determined automatically.

    two_tracers : boolean or None
        (Optional) boolean value sets 1- vs 2-tracer mode manually. If None (default), this mode is determined automatically.
    
    n_r_bins, n_mu_bins, n_jack, max_l : integer or None
        (Optional) integer value is used to set manually the number of radial bins, angular bins (not needed in some Legendre modes), jackknife regions (jackknife mode only) and maximum (even) ell (Legendre mode only) respectively. Each parameter which is None (default) is determined automatically.
    
    dry_run: boolean
        (Optional) If True, this will not run actual post-processing, only determine the filename and path (see below).

    Returns
    -------
    post_processing_results : dict[str, np.ndarray[float]]
        Post-processing results as a dictionary with string keys and Numpy array values. All this information is also saved in a ``Rescaled_Covariance_Matrices*.npz`` file in the ``out_dir`` (in ``file_root`` if the former is not provided).
        Selected common keys are: ``"full_theory_covariance"`` for the final covariance matrix and ``"shot_noise_rescaling"`` for the shot-noise rescaling value(s).
        For convenience, in the output dictionary only, ``"filename"`` contains the name of the file where the results were saved (which can be inconvenient to predict), and ``"path"`` contains its path (also obtainable by :func:`os.path.join`-ing ``out_dir`` with the filename)
    """
    # Set default output directory if not set
    if out_dir is None: out_dir = file_root

    # Determine mock post-processing
    mocks_precomputed = xi_sample_cov is not None
    mocks_from_samples = xi_11_samples is not None
    if mocks_precomputed and mocks_from_samples: raise ValueError("Provide either xi sample covariance or the xi samples, not both")
    if mocks_precomputed and load_sample_cov: raise ValueError("Either provide xi sample covariance or enable loading the sample covariance, not both")
    if mocks_from_samples and load_sample_cov: raise ValueError("Either provide xi samples or enable loading the sample covariance, not both")
    mocks_new = mocks_precomputed or mocks_from_samples

    jackknife_present = os.path.isdir(os.path.join(file_root, "xi_jack")) # whether the directory contains jackknife results. we might not use them depending on jackknife value

    if jackknife: # check for explicit incompatibilities
        if mocks_precomputed: raise ValueError("Provided xi sample covariance is not compatible with enabled jackknife")
        if mocks_from_samples: raise ValueError("Provided xi samples are not compatible with enabled jackknife")
        if load_sample_cov: raise ValueError("Loading xi sample covariance is not compatible with enabled jackknife")
    elif jackknife is None: # auto-determine jackknife
        if mocks_new or load_sample_cov: jackknife = False # jackknife must be disabled if (mock) sample covariance is set explicitly
        else: jackknife = jackknife_present

    # Auto-determine multi-tracer
    if two_tracers is None: two_tracers = os.path.isfile(os.path.join(file_root, f"xi/xi_22.dat"))

    ntracers = 2 if two_tracers else 1
    ncorr = ntracers * (ntracers + 1) // 2
    indices_corr = indices_corr_all[:ncorr]

    # Determine Legendre
    legendre_orig = len(glob(os.path.join(file_root, f"BinCorrectionFactor*"))) > 0
    legendre_mix = len(glob(os.path.join(file_root, "weights/mu_bin_legendre_factors_*.txt"))) > 0
    if legendre and not legendre_orig and not legendre_mix: warn("Legendre mode enabled explicitly, but characteristic Legendre files are not found")
    elif legendre is None: legendre = legendre_orig or legendre_mix

    print_function(f"Legendre: {legendre}")
    print_function(f"Jackknife: {jackknife}")
    print_function(f"Tuning to the provided (mock) sample covariance: {mocks_new}")
    print_function(f"Number of tracers: {ntracers}")

    if legendre and legendre_orig and jackknife: raise ValueError("Direct accumulation Legendre mode is not compatible with jackknives")

    # Determine number of radial, mu bins and/or jackknives automatically as needed
    binned_pair_names = glob("binned_pair_counts_n*_m*_j*_??.dat" if jackknife_present else "RR_counts_n*_m*_??.dat", root_dir = os.path.join(file_root, "weights"))
    if len(binned_pair_names) < ncorr: raise ValueError(f"Need {ncorr} pair counts, found {len(binned_pair_names)}")
    rstr = r'binned_pair_counts_n(?P<N_R_BINS>\d+)_m(?P<N_MU_BINS>\d+)_j(?P<N_JACK>\d+)_(?P<CORR_INDEX>\d+).dat' if jackknife_present else r'RR_counts_n(?P<N_R_BINS>\d+)_m(?P<N_MU_BINS>\d+)_(?P<CORR_INDEX>\d+).dat' # regex
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
    
    # Determine load_sample_cov automatically if needed
    mock_cov_basename = f"cov_sample_n{n_r_bins}_" + (f"l{max_l}" if legendre else f"m{n_mu_bins}") + ".txt"
    mock_cov_name = os.path.join(file_root, mock_cov_basename)
    if load_sample_cov:
        if not os.path.isfile(mock_cov_name): raise ValueError(f"The default (mock) sample covariance file, '{mock_cov_name}' does not exist and can not be loaded. Thus `load_sample_cov` can not be `True`.")
    elif load_sample_cov is None:
        if jackknife or mocks_new: load_sample_cov = False
        else: load_sample_cov = os.path.isfile(mock_cov_name)

    print_function(f"Tuning to the (mock) sample covariance loaded from the default file: {load_sample_cov}")

    mocks = mocks_new or load_sample_cov

    if not (jackknife or mocks):
        # cases when the shot-noise rescaling is not tuned - as it should be, or due to the lack of implementation
        print_function(f"Using {shot_noise_rescaling1=}" + two_tracers * f" and {shot_noise_rescaling2=}")

    if mocks_new:
        mock_cov_name = os.path.join(out_dir, mock_cov_basename) # in this case, the sample covariance will be written, and that should be into the output directory and not file_root; they can be the same if desired
        if os.path.isfile(mock_cov_name): warn(f"The default (mock) sample covariance file '{mock_cov_name}' exists and will be overwritten")
        elif not os.path.isdir(out_dir): # Need to make sure that the output directory exists. This is also checked in each post-processing functions, but only after writing the sample covariance file
            os.makedirs(os.path.abspath(out_dir)) # abspath is to exclude "../" for makedirs not to become "confused"
    # then write the sample covariance to file
    if mocks_precomputed:
        np.savetxt(mock_cov_name, xi_sample_cov)
    elif mocks_from_samples:
        if len(xi_11_samples) <= 1: raise ValueError("Need more than 1 sample in xi_11_samples to compute the sample covariance")
        if any(not np.allclose(xi_sample.edges[0], xi_11_samples[0].edges[0]) for xi_sample in xi_11_samples[1:]): raise ValueError(f"Found different separation/radial binning in xi_11_samples")
        if any(not np.allclose(xi_sample.edges[1], np.linspace(0, 1, n_mu_bins + 1)) for xi_sample in xi_11_samples): raise ValueError(f"Found angular/µ binning in xi_11_samples inconsistent with uniform and/or the inferred number of angular/µ bins")
        xi_samples_all = [xi_11_samples]
        if two_tracers:
            # check the required 22 samples
            if xi_22_samples is None: raise TypeError("xi_22_samples must be provided for multi-tracer")
            if len(xi_22_samples) != len(xi_11_samples): raise ValueError("xi_22_samples must contain the same number of samples as xi_11_samples")
            if any(not np.allclose(xi_sample.edges, xi_11_samples[0].edges) for xi_sample in xi_22_samples): raise ValueError(f"xi_22_samples must have the same binning as xi_11_samples")
            # check 12 samples which are not critical. if anything is wrong, substitute 11 samples as a placeholder
            if xi_12_samples is None:
                warn("xi_12_samples not provided. The shot-noise calibration should be fine, but do not use the sample covariance from the output folder directly.")
                xi_12_samples = xi_11_samples
            elif len(xi_12_samples) != len(xi_11_samples):
                warn("xi_12_samples must contain the same number of samples as xi_11_samples, will substitute them with a placeholder. The shot-noise calibration should be fine, but do not use the sample covariance from the output folder directly.")
                xi_12_samples = xi_11_samples
            elif any(not np.allclose(xi_sample.edges, xi_11_samples[0].edges) for xi_sample in xi_12_samples):
                warn(f"xi_12_samples must have the same binning as xi_11_samples, will substitute them with a placeholder. The shot-noise calibration should be fine, but do not use the sample covariance from the output folder directly.")
                xi_12_samples = xi_11_samples
            xi_samples_all += [xi_12_samples, xi_22_samples]
        if legendre:
            sample_cov_multipoles_from_pycorr_to_file(xi_samples_all, mock_cov_name, max_l=max_l)
        else:
            sample_cov_from_pycorr_to_file(xi_samples_all, mock_cov_name)

    if two_tracers:
        if legendre:
            if mocks:
                results = post_process_legendre_mocks_multi(mock_cov_name, file_root, n_r_bins, max_l, out_dir, skip_s_bins, skip_l, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
            elif jackknife:
                results = post_process_legendre_mix_jackknife_multi(*xi_jack_names, os.path.join(file_root, "weights"), file_root, n_mu_bins, max_l, out_dir, skip_s_bins, skip_l, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
            else:
                results = post_process_legendre_multi(file_root, n_r_bins, max_l, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_s_bins, skip_l, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        elif mocks:
            results = post_process_default_mocks_multi(mock_cov_name, file_root, n_r_bins, n_mu_bins, out_dir, skip_s_bins, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        elif jackknife:
            results = post_process_jackknife_multi(*xi_jack_names, os.path.join(file_root, "weights"), file_root, n_mu_bins, out_dir, skip_s_bins, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        else: # default
            results = post_process_default_multi(file_root, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_s_bins, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
    else:
        if legendre:
            if mocks:
                results = post_process_legendre_mocks(mock_cov_name, file_root, n_r_bins, max_l, out_dir, skip_s_bins, skip_l, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
            elif jackknife:
                results = post_process_legendre_mix_jackknife(xi_jack_names[0], os.path.join(file_root, "weights"), file_root, n_mu_bins, max_l, out_dir, skip_s_bins, skip_l, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
            else:
                results = post_process_legendre(file_root, n_r_bins, max_l, out_dir, shot_noise_rescaling1, skip_s_bins, skip_l, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        elif mocks:
            results = post_process_default_mocks(mock_cov_name, file_root, n_r_bins, n_mu_bins, out_dir, skip_s_bins, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        elif jackknife:
            results = post_process_jackknife(xi_jack_names[0], os.path.join(file_root, "weights"), file_root, n_mu_bins, out_dir, skip_s_bins, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)
        else: # default
            results = post_process_default(file_root, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, skip_s_bins, tracer = tracer, n_samples = n_samples, print_function = print_function, dry_run = dry_run)

    if extra_convergence_check and not dry_run:
        print_function("Performing an extra convergence check")
        results["extra_convergence_check_results"] = convergence_check_extra(results, print_function = print_function)

    return results
    