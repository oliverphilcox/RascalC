"""Implements an interface to estimate the covariance of 2-point correlation function."""

import numpy as np

def run_cov(mode: str, s_edges,
            nthread: int, N2: int, N3: int, N4: int, n_loops: int, loops_per_sample: int,
            out_dir: str, tmp_dir: str,
            randoms_positions1, randoms_weights1,
            pycorr_allcounts_11,
            xi_table_11, xi_table_12 = None, xi_table_22 = None,
            xi_cut_s: float = 250,
            pycorr_allcounts_12 = None, pycorr_allcounts_22 = None,
            no_data_galaxies1: float | None = None, no_data_galaxies2: float | None = None,
            randoms_samples1 = None,
            randoms_positions2 = None, randoms_weights2 = None, randoms_samples2 = None,
            n_mu_bins: int | None = None, max_l: int | None = None,
            boxsize: float | None = None):
    r"""
    Run the 2-point correlation function covariance integration.

    Parameters
    ----------
    mode : string
        Choice of binning setup, one of:

            - "default": compute covariance of the correlation function in s, µ bins. Only linear µ binning between 0 and 1 supported.
            - "legendre_mix": compute covariance of the correlation function Legendre multipoles in separation (s) bins projected from µ bins (only linear µ binning supported between 0 and 1). Works with jackknives, may be less efficient in periodic geometry.
            - "legendre_orig": compute covariance of the correlation function Legendre multipoles in separation (s) bins accumulated directly, without first doing µ-binned counts. Incompatible with jackknives.

    s_edges : sequence (list, array, tuple, etc) of floats
        Edges of the separation (s) bins.

    n_mu_bins : integer
        Number of µ bins (required in "default" ``mode``).
        Only linear µ binning between 0 and 1 supported.
    
    max_l : integer
        Max Legendre multipole index (required in both "legendre" ``mode``s).
        Has to be even.
    
    boxsize : None or float
        Periodic box side (one number only cubic supported so far).
        All the coordinates need to be between 0 and ``boxsize``.
        If None (default), assumed aperiodic.
    
    randoms_positions1 : array of floats, shape (3, N_randoms)
        Cartesian coordinates of random points for the first tracer.
    
    randoms_weights1 : array of floats of length N_randoms
        Weights of random points for the first tracer.
    
    randoms_samples1 : None or array of floats of length N_randoms
        (Optional) jackknife region numbers for random points for the first tracer.
        If given and not None, enables the jackknife functionality (tuning of shot-noise rescaling on jackknife correlation function estimates).
        The jackknife assignment must match the jackknife counts in ``pycorr_allcounts_11`` (and ``pycorr_allcounts_12`` in multi-tracer mode).

    randoms_positions2 : None or array of floats, shape (3, N_randoms2)
        (Optional) cartesian coordinates of random points for the second tracer.
        If given and not None, enables the multi-tracer functionality (full two-tracer covariance estimation).
    
    randoms_weights2 : None or array of floats of length N_randoms2
        Weights of random points for the second tracer (required for multi-tracer functionality).
    
    randoms_samples2 : None or array of floats of length N_randoms2
        Jackknife region numbers for the second tracer (required for multi-tracer + jackknife functionality, although this combination has not been used yet).
        The jackknife assignment must match the jackknife counts in ``pycorr_allcounts_12`` and ``pycorr_allcounts_22``.

    pycorr_allcounts_11 : ``pycorr`` TwoPointEstimator
        ``pycorr`` TwoPointEstimator with auto-counts for the first tracer.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples1``.

    pycorr_allcounts_12 : ``pycorr`` TwoPointEstimator
        (Optional) ``pycorr`` TwoPointEstimator with cross-counts between the two tracers.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples1`` and ``randoms_samples2``.

    pycorr_allcounts_22 : ``pycorr`` TwoPointEstimator
        (Optional) ``pycorr`` TwoPointEstimator with auto-counts for the second tracer.
        For jackknife functionality, must contain jackknife RR counts and correlation function. The jackknife assigment must match ``randoms_samples2``
    
    no_data_galaxies1 : None or float
        (Optional) number of first tracer data (not random!) points for the covariance rescaling.
        If not given or None, the code will attempt to obtain it from ``pycorr_allcounts_11``.
    
    no_data_galaxies2 : None or float
        (Optional) number of second tracer data (not random!) points for the covariance rescaling.
        If not given or None, the code will attempt to obtain it from ``pycorr_allcounts_22``.
    
    xi_table_11 : sequence (tuple or list) of 3 elements: (s_values, mu_values, xi_values)
        Table of first tracer auto-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.

    xi_table_12 : 
        Table of the two tracer's cross-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.

    xi_table_22 : 
        Table of second tracer auto-correlation function in separation (s) and µ bins.
        The code will use it for interpolation in the covariance matrix integrals.
    
    xi_cut_s : float
        (Optional) separation value beyond which the correlation function is assumed to be zero for the covariance matrix integrals.
        Between the maximum separation from ``xi_table``s and ``xi_cut_s``, the correlation function is extrapolated as :math:`\propto s^{-4}`.

    nthread : integer
        Number of hyperthreads to use.
        Can not utilize more threads than ``n_loops``.
    
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
        The runtime roughly scales as O(N_randoms * N2 * N3 * N4 * n_loops / nthread).

    loops_per_sample : integer
        Number of loops to merge into one output sample.
        Must divide ``max_loops``.
        Recommended to keep the number of samples = ``n_loops / loops_per_sample`` roughly between 10 and 30.

    out_dir : string
        Directory for important outputs.
        Moderate disk space required (up to a few hundred megabytes), but increases with covariance matrix size and number of samples (see above).

    tmp_dir : string
        Directory for temporary files.
        More disk space required - needs to store all the input arrays in the current implementation.
    """

    assert mode in ("default", "legendre_orig", "legendre_mix"), "Given mode not supported"

    # Set mode flags
    legendre_orig = (mode == "legendre_orig")
    legendre_mix = (mode == "legendre_mix")
    legendre = legendre_orig or legendre_mix
    jackknife = randoms_samples1 is not None

    assert not (legendre_orig and jackknife), "Original Legendre mode is not compatible with jackknife"

    if legendre:
        assert max_l is not None, "Max ell must be provided in Legendre mode"
        assert max_l % 2 == 0, "Only even Legendre multipoles supported"
    else:
        assert n_mu_bins is not None, "Number of µ bins for the covariance matrix must be provided in default mode"

    # Set some other flags
    periodic = boxsize is not None
    two_tracers = randoms_positions2 is not None

    if two_tracers: # check that everything is set accordingly
        assert randoms_weights2 is not None, "Second tracer weights must be provided in two-tracer mode"
        if jackknife: # although this case has not been used so far
            assert randoms_samples2 is not None, "Second tracer jackknife region numbers must be provided in two-tracer jackknife mode"
        assert pycorr_allcounts_12 is not None, "Cross-counts must be provided in two-tracer mode"
        assert pycorr_allcounts_22 is not None, "Second tracer auto-counts must be provided in two-tracer mode"
        assert xi_table_12 is not None, "Cross-correlation function must be provided in two-tracer mode"
        assert xi_table_22 is not None, "Second tracer auto-correlation function must be provided in two-tracer mode"
    
    assert n_loops % loops_per_sample == 0, "The sample collapsing factor must divide the number of loops"

    # Select the executable name
    exec_name = "bin/cov." + mode + "_jackknife" * jackknife + "_periodic" * periodic

    pass