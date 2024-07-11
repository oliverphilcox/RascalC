"""Implements an interface to estimate the covariance of 2-point correlation function."""

import pycorr
import numpy as np
import os
from datetime import datetime
from .pycorr_utils.utils import fix_bad_bins_pycorr, write_xi_file
from .write_binning_file import write_binning_file
from .pycorr_utils.jack import get_jack_xi_weights_counts_from_pycorr
from .pycorr_utils.counts import get_counts_from_pycorr
from .pycorr_utils.input_xi import get_input_xi_from_pycorr
from .mu_bin_legendre_factors import write_mu_bin_legendre_factors
from .correction_function import compute_correction_function, compute_correction_function_multi
from .convergence_check_extra import convergence_check_extra
from .utils import rmdir_if_exists_and_empty


suffixes_tracer_all = ("", "2") # all supported tracer suffixes
indices_corr_all = ("11", "12", "22") # all supported 2PCF indices
suffixes_corr_all = ("", "12", "2") # all supported 2PCF suffixes
tracer1_corr = (0, 0, 1)
tracer2_corr = (0, 1, 1)


def run_cov(mode: str,
            nthread: int, N2: int, N3: int, N4: int, n_loops: int, loops_per_sample: int,
            out_dir: str, tmp_dir: str,
            randoms_positions1: np.ndarray[float], randoms_weights1: np.ndarray[float],
            pycorr_allcounts_11: pycorr.twopoint_estimator.BaseTwoPointEstimator,
            xi_table_11: pycorr.twopoint_estimator.BaseTwoPointEstimator | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]] | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]],
            position_type: str = "pos",
            xi_table_12: None | pycorr.twopoint_estimator.BaseTwoPointEstimator | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]] | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] = None,
            xi_table_22: None | pycorr.twopoint_estimator.BaseTwoPointEstimator | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]] | tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] = None,
            xi_cut_s: float = 250, xi_refinement_iterations: int = 10,
            pycorr_allcounts_12: pycorr.twopoint_estimator.BaseTwoPointEstimator | None = None, pycorr_allcounts_22: pycorr.twopoint_estimator.BaseTwoPointEstimator | None = None,
            normalize_wcounts: bool = True,
            no_data_galaxies1: float | None = None, no_data_galaxies2: float | None = None,
            randoms_samples1: np.ndarray[int] | None = None,
            randoms_positions2: np.ndarray[float] | None = None, randoms_weights2: np.ndarray[float] | None = None, randoms_samples2: np.ndarray[int] | None = None,
            max_l: int | None = None,
            boxsize: float | None = None,
            skip_s_bins: int = 0, skip_l: int = 0,
            shot_noise_rescaling1: float = 1, shot_noise_rescaling2: float = 1,
            sampling_grid_size: int = 301, coordinate_scaling: float = 1, seed: int | None = None,
            verbose: bool = False) -> dict[str, np.ndarray[float]]:
    r"""
    Run the 2-point correlation function covariance integration.

    Parameters
    ----------
    mode : string
        Choice of binning setup, one of:

            - "s_mu": compute covariance of the correlation function in s, µ bins. Only linear µ binning between 0 and 1 supported.
            - "legendre_projected": compute covariance of the correlation function Legendre multipoles in separation (s) bins projected from µ bins (only linear µ binning supported between 0 and 1). Procedure matches ``pycorr``. Works with jackknives, may be less efficient in periodic geometry.
            - "legendre_accumulated": compute covariance of the correlation function Legendre multipoles in separation (s) bins accumulated directly, without first doing µ-binned counts. Incompatible with jackknives.
    
    max_l : integer
        Max Legendre multipole index (required in both "legendre" ``mode``\s).
        Has to be even.
    
    boxsize : None or float
        Periodic box side (one number only cubic supported so far).
        All the coordinates need to be between 0 and ``boxsize``.
        If None (default), assumed aperiodic.

    position_type : string, default="pos"
        Type of input positions, one of:

            - "rdd": RA, Dec in degrees, distance
            - "xyz": Cartesian positions, shape (3, N)
            - "pos": Cartesian positions, shape (N, 3).
    
    randoms_positions1 : array of floats, shaped according to `position_type`
        Cartesian coordinates of random points for the first tracer.
    
    randoms_weights1 : array of floats of length N_randoms
        Weights of random points for the first tracer.
    
    randoms_samples1 : None or array of floats of length N_randoms
        (Optional) jackknife region numbers for random points for the first tracer.
        If given and not None, enables the jackknife functionality (tuning of shot-noise rescaling on jackknife correlation function estimates).
        The jackknife assignment must match the jackknife counts in ``pycorr_allcounts_11`` (and ``pycorr_allcounts_12`` in multi-tracer mode).

    randoms_positions2 : None or array of floats, shaped according to `position_type`
        (Optional) cartesian coordinates of random points for the second tracer.
        If given and not None, enables the multi-tracer functionality (full two-tracer covariance estimation).
    
    randoms_weights2 : None or array of floats of length N_randoms2
        Weights of random points for the second tracer (required for multi-tracer functionality).
    
    randoms_samples2 : None or array of floats of length N_randoms2
        Jackknife region numbers for the second tracer (required for multi-tracer + jackknife functionality, although this combination has not been used yet).
        The jackknife assignment must match the jackknife counts in ``pycorr_allcounts_12`` and ``pycorr_allcounts_22``.

    pycorr_allcounts_11 : ``pycorr.TwoPointEstimator``
        ``pycorr.TwoPointEstimator`` with auto-counts for the first tracer.
        Must be rebinned and/or cut to the separation (s) bins desired for the covariance.
        Note that more bins result in slower convergence. A typical configuration has been 4 Mpc/h wide bins between 20 and 200 Mpc/h.
        The counts will be wrapped to positive µ, so if the µ range in them is from -1 to 1, the number of µ bins must be even.
        Providing unwrapped counts (µ from -1 to 1) is preferrable, because some issues can be fixed by assuming symmetry.
        In "s_mu" ``mode``, the covariance will be done for the given number of µ bins (after wrapping).
        In "legendre_projected" ``mode``, it will be assumed that Legendre multipoles are projected from the same number of µ bins as present in these counts. One might consider rebinning more coarsely for faster performance but less guaranteed accuracy (neither effect has been tested yet).
        In "legendre_accumulated" ``mode``, all the present µ bins (after wrapping) will be used to fit the survey correction functions.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples1``.

    pycorr_allcounts_12 : ``pycorr.TwoPointEstimator``
        (Optional) ``pycorr.TwoPointEstimator`` with cross-counts between the two tracers (order does not matter, because they will be wrapped).
        Must have the same bin configuration as ``pycorr_allcounts_11``.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples1`` and ``randoms_samples2``.

    pycorr_allcounts_22 : ``pycorr.TwoPointEstimator``
        (Optional) ``pycorr.TwoPointEstimator`` with auto-counts for the second tracer.
        Must have the same bin configuration as ``pycorr_allcounts_11``.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples2``

    normalize_wcounts : boolean
        (Optional) whether to normalize the weights and weighted counts.
        If False, the provided RR counts must match what can be obtained from given randoms, otherwise the covariance matrix will be off by a constant factor.
        Example: if counts were computed with ``n_randoms`` roughly similar random chunks and only one is provided to RascalC here, the counts should be divided by ``n_random`` where ``s > split_above`` and by ``n_random ** 2`` where ``s < split_above``.
        If True (default), the weights will be normalized so that their sum is 1 and the counts will be normalized by their ``wnorm``, which gives a match with default ``pycorr`` normalization settings.
    
    no_data_galaxies1 : None or float
        (Optional) number of first tracer data (not random!) points for the covariance rescaling.
        If not given or None, the code will attempt to obtain it from ``pycorr_allcounts_11``.
    
    no_data_galaxies2 : None or float
        (Optional) number of second tracer data (not random!) points for the covariance rescaling.
        If not given or None, the code will attempt to obtain it from ``pycorr_allcounts_22``.
    
    xi_table_11 : ``pycorr.TwoPointEstimator``, or sequence (tuple or list) of 3 elements: (s_values, mu_values, xi_values), or sequence (tuple or list) of 4 elements: (s_values, mu_values, xi_values, s_edges)
        Table of first tracer auto-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.
        Important: if the given correlation function is an average in s, µ bins, the separation bin edges need to be provided (and the µ bins are assumed to be linear) for rescaling procedure which ensures that the interpolation results averaged over s, µ bins returns the given correlation function. In case of ``pycorr.TwoPointEstimator``, the edges will be recovered automatically. Unwrapped estimators (µ from -1 to 1) are preferred, because symmetry allows to fix some issues.
        In the sequence format:

            - s_values must be a 1D array of reference separation (s) values for the table, of length N;
            - mu_values must be a 1D array of reference µ values (covering the range from 0 to 1) for the table, of length M;
            - xi_values must be an array of correlation function values at those s, µ values of shape (N, M);
            - s_edges, if given, must be a 1D array of separation bin edges of length N+1. The bins must come close to zero separation (say start at ``s <= 0.01``).
        
        The sequence containing 3 elements should be used for theoretical models evaluated at a grid of s, mu values.
        The 4-element format should be used for bin-averaged estimates.

    xi_table_12 : None or the same format as ``xi_table_11``
        Table of the two tracer's cross-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.
        Required for multi-tracer functionality.

    xi_table_22 : None or the same format as ``xi_table_11``
        Table of second tracer auto-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.
        Required for multi-tracer functionality.
    
    xi_cut_s : float
        (Optional) separation value beyond which the correlation function is assumed to be zero for the covariance matrix integrals. Default: 250.
        Between the maximum separation from ``xi_table``s and ``xi_cut_s``, the correlation function is extrapolated as :math:`\propto s^{-4}`.
    
    xi_refinement_iterations : int
        (Optional) number of iterations in the correlation function refinement procedure for interpolation inside the code, ensuring that the bin-averaged interpolated values match the binned correlation function estimates. Default: 10.
        Important: the refinement procedure is disabled completely regardless of this setting if the ``xi_table``\s are sequences of 3 elements, since they are presumed to be a theoretical model evaluated at a grid of s, mu values and not bin averages.

    nthread : integer
        Number of hyperthreads to use.
        Can not utilize more threads than ``n_loops``.
        IMPORTANT: AVOID multi-threading in the Python process calling this function (e.g. at NERSC this would mean not setting `OMP_*` and other `*_THREADS` environment variables; the code should be able to set them by itself). Otherwise the code may run effectively single-threaded. If you need other multi-threaded calculations, run them separately or spawn sub-processes.
    
    N2 : integer
        Number of secondary points to sample per each primary random point.
        Setting too low (below 5) is not recommended.
    
    N3 : integer
        Number of tertiary points to sample per each secondary point.
        Setting too low (below 5) is not recommended.
    
    N4 : integer
        Number of quaternary points to sample per each tertiary point.
        Setting too low (below 5) is not recommended.

    n_loops : integer
        Number of integration loops.
        For optimal balancing and minimal idle time, should be a few times (at least twice) ``nthread`` and exactly divisible by it.
        The runtime roughly scales as the number of quads per the number of threads, O(N_randoms * N2 * N3 * N4 * n_loops / nthread).
        For reference, on NERSC Perlmutter CPU node the code processed about 5 millions (5e6) quads per second per (hyper)thread (a node has 256 of them) as of October 2023. (In Legendre projected mode, which is probably the slowest, with N2=5, N3=10, N4=20.)
        In single-tracer mode, the number of quads is ``N_randoms * N2 * N3 * N4 * n_loops``.
        In two-tracer mode, the number of quads is ``(5 * N_randoms1 + 2 * N_randoms2) * N2 * N3 * N4 * n_loops``.

    loops_per_sample : integer
        Number of loops to merge into one output sample.
        Must divide ``max_loops``.
        Recommended to keep the number of samples = ``n_loops / loops_per_sample`` roughly between 10 and 30.

    out_dir : string
        Directory for important outputs.
        Moderate disk space required (up to a few hundred megabytes), but increases with covariance matrix size and number of samples (see above).
        Avoid ".." in this path because it makes os.makedirs() "confused".

    tmp_dir : string
        Directory for temporary files. Contents can be deleted after the code has run, but this will not be done automatically.
        More disk space required - needs to store all the input arrays in the current implementation.
        Avoid ".." in this path because it makes os.makedirs() "confused".

    skip_s_bins : int
        (Optional) number of lowest separations bins to skip at the post-processing stage. Those tend to converge worse and probably will not be precise due to the limitations of the formalism. Default 0 (no skipping).

    skip_l : int
        (Only for the Legendre modes; optional) number of highest (even) multipoles to skip at the post-processing stage. Those tend to converge worse. Default 0 (no skipping).

    shot_noise_rescaling1 : float
        (Optional) shot-noise rescaling value for the first tracer if known beforehand. Default 1 (no rescaling).
        Will be ignored in jackknife mode - then shot-noise rescaling is optimized on the auto-covariance.

    shot_noise_rescaling2 : float
        (Optional) shot-noise rescaling value for the second tracer if known beforehand. Default 1 (no rescaling).
        Will be ignored in jackknife mode - then shot-noise rescaling is optimized on the auto-covariance.

    sampling_grid_size : int
        (Optional) first guess for the sampling grid size.
        The code should be able to find a suitable number automatically.

    coordinate_scaling : float
        (Optional) scaling factor for all the Cartesian coordinates. Default 1 (no rescaling).
        This option is supported by the C++ code, but its use cases are unclear.

    seed : integer or None
        (Optional) If given as an integer, allows to reproduce the results with the same settings, except the number of threads.
        Individual subsamples may differ because they are accumulated/written in order of loop completion which may depend on external factors at runtime, but the final integrals should be the same.
        If None (default), the initialization will be random.

    verbose : bool
        (Optional) report each 5% of each loop's progress by printing. Default False (off).
        This can be a lot of output, only use when the number of loops is small.

    Returns
    -------
    post_processing_results : dict[str, np.ndarray[float]]
        Post-processing results as a dictionary with string keys and Numpy array values. All this information is also saved in a Rescaled_Covariance_Matrices*.npz file in the output directory.
        Selected common keys are: "full_theory_covariance" for the final covariance matrix and "shot_noise_rescaling" for the shot-noise rescaling value(s).
        There will also be a Raw_Covariance_Matices*.npz file in the output directory (as long as the C++ code has run without errors), which can be post-processed separately in a different way.
    """

    if mode not in ("s_mu", "legendre_accumulated", "legendre_projected"): raise ValueError("Given mode not supported")

    # Set mode flags
    legendre_orig = (mode == "legendre_accumulated")
    legendre_mix = (mode == "legendre_projected")
    legendre = legendre_orig or legendre_mix
    jackknife = randoms_samples1 is not None

    if (legendre_orig and jackknife): raise ValueError("Direct accumulation Legendre mode is not compatible with jackknives")

    if legendre:
        if max_l is None: raise TypeError("Max ell must be provided in Legendre mode")
        if max_l % 2 != 0: raise ValueError("Only even Legendre multipoles supported")
    else:
        if n_mu_bins is None: raise TypeError("Number of µ bins for the covariance matrix must be provided in s_mu mode")

    # Set some other flags
    periodic = bool(boxsize) # False for None (default) and 0
    two_tracers = randoms_positions2 is not None

    if two_tracers: # check that everything is set accordingly
        if randoms_weights2 is None: raise TypeError("Second tracer weights must be provided in two-tracer mode")
        if jackknife: # although this case has not been used so far
            if randoms_samples2 is None: raise TypeError("Second tracer jackknife region numbers must be provided in two-tracer jackknife mode")
        if pycorr_allcounts_12 is None: raise TypeError("Cross-counts must be provided in two-tracer mode")
        if pycorr_allcounts_22 is None: raise TypeError("Second tracer auto-counts must be provided in two-tracer mode")
        if xi_table_12 is None: raise TypeError("Cross-correlation function must be provided in two-tracer mode")
        if xi_table_22 is None: raise TypeError("Second tracer auto-correlation function must be provided in two-tracer mode")
    
    if n_loops % loops_per_sample != 0: raise ValueError("The sample collapsing factor must divide the number of loops")

    ntracers = 2 if two_tracers else 1
    ncorr = ntracers * (ntracers + 1) // 2
    suffixes_tracer = suffixes_tracer_all[:ntracers]
    indices_corr = indices_corr_all[:ncorr]
    suffixes_corr = suffixes_corr_all[:ncorr]

    s_edges = pycorr_allcounts_11.edges[0]
    n_r_bins = len(s_edges) - 1
    n_mu_bins = pycorr_allcounts_11.shape[1]
    if pycorr_allcounts_11.edges[1][0] < 0: n_mu_bins //= 2 # will be wrapped
    if jackknife:
        jack_region_numbers = np.unique(randoms_samples1)
        njack = len(jack_region_numbers)
        if two_tracers:
            if not np.array_equal(jack_region_numbers, np.unique(randoms_samples2)): # comparison is good because unique results are sorted
                raise ValueError("The sets of jackkknife labels of the two tracers must be the same")

    # set the technical filenames
    input_filenames = [os.path.join(tmp_dir, str(t) + ".txt") for t in range(ntracers)]
    cornames = [os.path.join(out_dir, f"xi/xi_{index}.dat") for index in indices_corr]
    binned_pair_names = [os.path.join(out_dir, "weights/" + ("binned_pair" if jackknife else "RR") + f"_counts_n{n_r_bins}_m{n_mu_bins}" + (f"_j{njack}" if jackknife else "") + f"_{index}.dat") for index in indices_corr]
    if jackknife:
        jackknife_weights_names = [os.path.join(out_dir, f"weights/jackknife_weights_n{n_r_bins}_m{n_mu_bins}_j{njack}_{index}.dat") for index in indices_corr]
        xi_jack_names = [os.path.join(out_dir, f"xi_jack/xi_jack_n{n_r_bins}_m{n_mu_bins}_j{njack}_{index}.dat") for index in indices_corr]
        jackknife_pairs_names = [os.path.join(out_dir, f"weights/jackknife_pair_counts_n{n_r_bins}_m{n_mu_bins}_j{njack}_{index}.dat") for index in indices_corr]
    if legendre_orig:
        phi_names = [os.path.join(out_dir, f"BinCorrectionFactor_n{n_r_bins}_" + ("periodic" if periodic else f'm{n_mu_bins}') + f"_{index}.txt") for index in indices_corr]
    
    # make sure the dirs exist
    os.makedirs(out_dir, exist_ok = True)
    os.makedirs(tmp_dir, exist_ok = True)
    os.makedirs(os.path.join(out_dir, "xi"), exist_ok = True)
    os.makedirs(os.path.join(out_dir, "weights"), exist_ok = True)
    if jackknife: os.makedirs(os.path.join(out_dir, "xi_jack"), exist_ok = True)
    
    # Create a log file in output directory
    logfilename = "log.txt"
    logfile = os.path.join(out_dir, logfilename)

    def print_log(s: object) -> None:
        os.system(f"echo \"{s}\" >> {logfile}")

    def print_and_log(s: object) -> None:
        print(s)
        print_log(s)
    
    print_and_log(f"Mode: {mode}")
    print_and_log(f"Periodic box: {periodic}")
    if periodic: print_and_log(f"Box side: {boxsize}")
    print_and_log(f"Jackknife: {jackknife}")
    print_and_log(f"Number of tracers: {1 + two_tracers}")
    print_and_log(f"Normalizing weights and weighted counts: {normalize_wcounts}")
    print_and_log(datetime.now())

    counts_factor = None if normalize_wcounts else 1
    ndata = (no_data_galaxies1, no_data_galaxies2)[:ntracers]

    # convert counts and jackknife xi if needed; loading ndata too whenever not given
    pycorr_allcounts_all = (pycorr_allcounts_11, pycorr_allcounts_12, pycorr_allcounts_22)
    for c, pycorr_allcounts in enumerate(pycorr_allcounts_all[:ncorr]):
        if pycorr_allcounts.edges[1][0] < 0: # try to fix and wrap
            pycorr_allcounts = fix_bad_bins_pycorr(pycorr_allcounts)
            print_and_log(f"Wrapping pycorr_allcounts_{indices_corr[c]} to µ>=0")
            pycorr_allcounts = pycorr_allcounts.wrap()
        if not np.allclose(pycorr_allcounts.edges[1], np.linspace(0, 1, n_mu_bins + 1)): raise ValueError(f"pycorr_allcounts_{indices_corr[c]} µ binning is not consistent with linear between 0 and 1 (after wrapping)")
        RR_counts = get_counts_from_pycorr(pycorr_allcounts, counts_factor)
        np.savetxt(binned_pair_names[c], RR_counts.reshape(-1, 1)) # the file needs to have 1 column
        if jackknife:
            xi_jack, jack_weights, jack_RR_counts = get_jack_xi_weights_counts_from_pycorr(pycorr_allcounts, counts_factor)
            if not np.allclose(np.sum(jack_RR_counts, axis=0), RR_counts.ravel()): raise ValueError("Total counts mismatch")
            ## Write to files using numpy functions
            write_xi_file(xi_jack_names[c], pycorr_allcounts.sepavg(axis = 0), pycorr_allcounts.sepavg(axis = 1), xi_jack)
            jack_numbers = pycorr_allcounts.realizations # column of jackknife numbers, may be useless but needed for format compatibility
            if not np.array_equal(np.array(jack_region_numbers, dtype = int), np.sort(np.array(jack_numbers, dtype = int))): raise ValueError(f"The code requires integer jackknife numbers consistent between the randoms_samples1 and pycorr_allcounts_{indices_corr[c]}") # jack_region_numbers are the unique jackknife labels from the data (already sorted)
            np.savetxt(jackknife_weights_names[c], np.column_stack((jack_numbers, jack_weights)))
            np.savetxt(jackknife_pairs_names[c], np.column_stack((jack_numbers, jack_RR_counts))) # not really needed for the C++ code or processing but let it be
        # fill ndata if not given
        tracer1 = tracer1_corr[c]
        if not ndata[tracer1]:
            ndata[tracer1] = pycorr_allcounts.D1D2.size1

    if any(not tracer_ndata for tracer_ndata in ndata): raise ValueError("Not given and not recovered all the necessary normalization factors (no_data_galaxies1/2)")
    print_and_log(f"Number(s) of data galaxies: {ndata}")
    
    # write the xi file(s); need to set the 2PCF binning (even if only technical) and decide whether to rescale the 2PCF in the C++ code
    all_xi = (xi_table_11, xi_table_12, xi_table_22)
    xi_s_edges = None
    xi_n_mu_bins = None
    refine_xi = False
    for c, xi in enumerate(all_xi[:ncorr]):
        if c > 0:
            if type(xi) != type(all_xi[0]): raise TypeError(f"xi_table_{indices_corr[c]} must have the same type as xi_table_11")
            if len(xi) != len(all_xi[0]): raise ValueError(f"xi_table_{indices_corr[c]} must have the same structure as xi_table_11")
        if isinstance(xi, pycorr.twopoint_estimator.BaseTwoPointEstimator):
            refine_xi = True
            if xi.edges[1][0] < 0:
                xi = fix_bad_bins_pycorr(xi)
                print_and_log(f"Wrapping xi_table_{indices_corr[c]} to µ>=0")
                xi = xi.wrap()
            if c == 0:
                xi_n_mu_bins = xi.shape[1]
                xi_s_edges = xi.edges[0]
            elif not np.allclose(xi_s_edges, xi.edges[0]): raise ValueError("Different binning for different correlation functions not supported")
            if not np.allclose(xi.edges[1], np.linspace(0, 1, xi_n_mu_bins + 1)): raise ValueError(f"xi_table_{indices_corr[c]} µ binning is not consistent with linear between 0 and 1 (after wrapping)")
            write_xi_file(cornames[c], xi.sepavg(axis = 0), xi.sepavg(axis = 1), get_input_xi_from_pycorr(xi))
        elif isinstance(xi, tuple) or isinstance(xi, list):
            if len(xi) == 4: # the last element is the edges
                refine_xi = True
                if c == 0: xi_s_edges = xi[-1]
                elif not np.allclose(xi_s_edges, xi[-1]): raise ValueError("Different binning for different correlation functions not supported")
                xi = xi[:-1]
            if len(xi) != 3: raise ValueError(f"xi_table {indices_corr[c]} must have 3 or 4 elements if a tuple/list")
            r_vals, mu_vals, xi_vals = xi
            if len(xi_vals) != len(r_vals): raise ValueError(f"xi_values {indices_corr[c]} must have the same number of rows as r_values")
            if len(xi_vals[0]) != len(mu_vals): raise ValueError(f"xi_values {indices_corr[c]} must have the same number of columns as mu_values")
            if c == 0:
                xi_n_mu_bins = len(mu_vals)
                if not refine_xi:
                    xi_s_edges = (r_vals[:-1] + r_vals[1:]) / 2 # middle values as midpoints of r_vals to be safe
                    xi_s_edges = [1e-4] + xi_s_edges + [2 * r_vals[-1] - xi_s_edges[-1]] # set the lowest edge near 0 and the highest beyond the last point of r_vals
            write_xi_file(cornames[c], r_vals, mu_vals, xi_vals)
        else: raise TypeError(f"Xi table {indices_corr[c]} must be either a pycorr.TwoPointEstimator or a tuple/list")
    xi_refinement_iterations *= refine_xi # True is 1; False is 0 => 0 iterations => no refinement
    
    # write the randoms file(s)
    randoms_positions = [randoms_positions1, randoms_positions2]
    randoms_weights = [randoms_weights1, randoms_weights2]
    randoms_samples = (randoms_samples1, randoms_samples2)
    for t, input_filename in enumerate(input_filenames):
        randoms_properties = pycorr.twopoint_counter._format_positions(randoms_positions[t], mode = "smu", position_type = position_type, dtype = np.float64) # list of x, y, z coordinate arrays; weights (and jackknife region numbers if any) will be appended
        if legendre_orig: randoms_positions[t] = np.array(randoms_properties) # save the formatted positions as an array for correction function computation
        nrandoms = len(randoms_properties[0])
        if randoms_weights[t].ndim != 1: raise ValueError(f"Weights of randoms {t+1} not contained in a 1D array")
        if len(randoms_weights[t]) != nrandoms: raise ValueError(f"Number of weights for randoms {t+1} mismatches the number of positions")
        if normalize_wcounts: randoms_weights[t] /= np.sum(randoms_weights[t])
        randoms_properties.append(randoms_weights[t])
        if jackknife:
            if randoms_samples[t].ndim != 1: raise ValueError(f"Weights of sample labels {t+1} not contained in a 1D array")
            if len(randoms_samples[t]) != nrandoms: raise ValueError(f"Number of sample labels for randoms {t+1} mismatches the number of positions")
            randoms_properties.append(randoms_samples[t])
        np.savetxt(input_filename, np.column_stack(randoms_properties))
        randoms_properties = None

    # write the binning files
    binfile = os.path.join(out_dir, "radial_binning_cov.csv")
    write_binning_file(binfile, s_edges)
    binfile_cf = os.path.join(out_dir, "radial_binning_corr.csv")
    write_binning_file(binfile_cf, xi_s_edges)

    # Select the executable name
    exec_name = "bin/cov." + mode + "_jackknife" * jackknife + "_periodic" * periodic + "_verbose" * verbose
    # the above must be true relative to the script location
    # below we should make it absolute, i.e. right regardless of the working directory
    exec_path = os.path.join(os.path.realpath(os.path.dirname(__file__)), exec_name)

    # form the command line
    command = "env OMP_PROC_BIND=spread OMP_PLACES=threads " # set OMP environment variables, should not be set before
    command += f"{exec_path} -output {out_dir}/ -nside {sampling_grid_size} -rescale {coordinate_scaling} -nthread {nthread} -maxloops {n_loops} -loopspersample {loops_per_sample} -N2 {N2} -N3 {N3} -N4 {N4} -xicut {xi_cut_s} -binfile {binfile} -binfile_cf {binfile_cf} -mbin_cf {xi_n_mu_bins} -cf_loops {xi_refinement_iterations}" # here are universally acceptable parameters
    command += "".join([f" -in{suffixes_tracer[t]} {input_filenames[t]}" for t in range(ntracers)]) # provide all the random filenames
    command += "".join([f" -norm{suffixes_tracer[t]} {ndata[t]}" for t in range(ntracers)]) # provide all ndata for normalization
    command += "".join([f" -cor{suffixes_corr[c]} {cornames[c]}" for c in range(ncorr)]) # provide all correlation functions
    if legendre: # only provide max multipole l for now
        command += f" -max_l {max_l}"
    if legendre_mix: # generate and provide factors filename
        mu_bin_legendre_file = write_mu_bin_legendre_factors(n_mu_bins, max_l, os.path.join(out_dir, "weights"))
        command += f" -mu_bin_legendre_file {mu_bin_legendre_file}"
    if not legendre_orig: # provide binned pair counts files and number of mu bin
        command += "".join([f" -RRbin{suffixes_corr[c]} {binned_pair_names[c]}" for c in range(ncorr)]) + f" -mbin {n_mu_bins}"
    if periodic: # append periodic flag and box size
        command += f" -perbox -boxsize {boxsize}"
    if jackknife: # provide jackknife weight files for all correlations
        command += "".join([f" -jackknife{suffixes_corr[c]} {jackknife_weights_names[c]}" for c in range(ncorr)])

    # compute the correction function if original Legendre
    if legendre_orig: # need correction function
        print_and_log(datetime.now())
        print_and_log(f"Computing the correction function")
        if ntracers == 1:
            compute_correction_function(randoms_positions[0], randoms_weights[0], binfile, out_dir, periodic, binned_pair_names[0], print_and_log)
        elif ntracers == 2:
            compute_correction_function_multi(randoms_positions[0], randoms_weights[0], randoms_positions[1], randoms_weights[1], binfile, out_dir, periodic, *binned_pair_names, print_function = print_and_log)
        command += "".join([f" -phi_file{suffixes_corr[c]} {phi_names[c]}" for c in range(ncorr)])

    # deal with the seed
    if seed is not None: # need to pass to the C++ code and make sure it can be received properly
        if not isinstance(seed, int): raise TypeError("Seed must be int or None")
        seed &= 2**32 - 1 # this bitwise AND truncates the seed into a 32-bit unsigned (positive) integer (definitely a subset of unsigned long)
        command += f" -seed {seed}"
    
    # run the main code
    print_and_log(datetime.now())
    print_and_log(f"Launching the C++ code with command: {command}")
    status = os.system(f"bash -c 'set -o pipefail; stdbuf -oL -eL {command} 2>&1 | tee -a {logfile}'")
    # tee prints what it gets to stdout AND saves to file
    # stdbuf -oL -eL should solve the output delays due to buffering without hurting the performance too much
    # without pipefail, the exit_code would be of tee, not reflecting main command failures
    # feed the command to bash because on Ubuntu it was executed in sh (dash) where pipefail is not supported

    # clean up
    for input_filename in input_filenames: os.remove(input_filename) # delete the larger (temporary) input files
    rmdir_if_exists_and_empty(tmp_dir) # safely remove the temporary directory

    # check the run status
    exit_code = os.waitstatus_to_exitcode(status) # assumes we are in Unix-based OS; on Windows status is the exit code
    if exit_code: raise RuntimeError(f"The C++ code terminated with an error: exit code {exit_code}")
    print_and_log("The C++ code finished succesfully")

    # post-processing
    print_and_log(datetime.now())
    print_and_log("Starting post-processing")
    if two_tracers:
        if legendre:
            from .post_process import post_process_legendre_multi
            results = post_process_legendre_multi(out_dir, n_r_bins, max_l, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_s_bins, skip_l, print_function = print_and_log)
        elif jackknife:
            from .post_process import post_process_jackknife_multi
            results = post_process_jackknife_multi(*xi_jack_names, os.path.dirname(jackknife_weights_names[0]), out_dir, n_mu_bins, out_dir, skip_s_bins, print_function = print_and_log)
        else: # default
            from .post_process import post_process_default_multi
            results = post_process_default_multi(out_dir, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, shot_noise_rescaling2, skip_s_bins, print_function = print_and_log)
    else:
        if legendre:
            if jackknife:
                from .post_process import post_process_legendre_mix_jackknife
                results = post_process_legendre_mix_jackknife(xi_jack_names[0], os.path.dirname(jackknife_weights_names[0]), out_dir, n_mu_bins, max_l, out_dir, skip_s_bins, skip_l, print_function = print_and_log)
            else:
                from .post_process import post_process_legendre
                results = post_process_legendre(out_dir, n_r_bins, max_l, out_dir, shot_noise_rescaling1, skip_s_bins, skip_l, print_function = print_and_log)
        elif jackknife:
            from .post_process import post_process_jackknife
            results = post_process_jackknife(xi_jack_names[0], os.path.dirname(jackknife_weights_names[0]), out_dir, n_mu_bins, out_dir, skip_s_bins, print_function = print_and_log)
        else: # default
            from .post_process import post_process_default
            results = post_process_default(out_dir, n_r_bins, n_mu_bins, out_dir, shot_noise_rescaling1, skip_s_bins, print_function = print_and_log)

    print_and_log("Finished post-processing")
    print_and_log(datetime.now())

    print_and_log("Performing an extra convergence check")
    convergence_check_extra(results, print_function = print_and_log)

    print_and_log("Finished.")
    print_and_log(datetime.now())
    return results