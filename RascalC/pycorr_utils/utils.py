import pycorr
import numpy as np
from warnings import warn
from ..utils import format_skip_r_bins
from ..xi.utils import write_xi_file # for convenience of use in other modules


def fix_bad_bins_counts(counts: pycorr.twopoint_counter.BaseTwoPointCounter) -> pycorr.twopoint_counter.BaseTwoPointCounter:
    """
    Takes a counts object and fixes bins with negative wcounts by overwriting their content by reflection.
    Only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32.
    """
    mu_edges = counts.edges[1]
    if not np.allclose(mu_edges, -mu_edges[::-1]):
        raise ValueError(f'input counts cannot be fixed by symmetry as 2nd dimension edges are not symmetric: {mu_edges} != {-mu_edges[::-1]}')
    bad_bins_mask = counts.wcounts < 0
    for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
        warn(f"Negative {counts.name}.wcounts ({counts.wcounts[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.wcounts[s_bin, -1-mu_bin]:.2e})")
        counts.wcounts[s_bin, mu_bin] = counts.wcounts[s_bin, -1-mu_bin]
    if np.ndim(counts.wnorm) == 2: # can't check and fix wnorm if it is a scalar (i.e., the same for all bins) or a 1D array. previously, in most practical cases, wnorm was a 2D array because counts below 20 Mpc/h were from concatenated randoms and from split/disjoint randoms above 20 Mpc/h (in principle, in such cases it could be a 1D array barring the potential array broadcasting issues)
        bad_bins_mask = counts.wnorm < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {counts.name}.wnorm ({counts.wnorm[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.wnorm[s_bin, -1-mu_bin]:.2e})")
            counts.wnorm[s_bin, mu_bin] = counts.wnorm[s_bin, -1-mu_bin]
    return counts


def fix_bad_bins_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator) -> pycorr.twopoint_estimator.BaseTwoPointEstimator:
    """
    Fixes bins with negative wcounts by overwriting their content by reflection.
    Only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32.
    """
    cls = xi_estimator.__class__
    kw = {}
    for name in xi_estimator.count_names:
        counts : pycorr.twopoint_counter.BaseTwoPointCounter = getattr(xi_estimator, name)
        counts = fix_bad_bins_counts(counts)
        if isinstance(counts, pycorr.twopoint_jackknife.JackknifeTwoPointCounter): # need to fix the counts in all realizations of the jackknife counter
            for i in xi_estimator.realizations:
                counts.auto[i] = fix_bad_bins_counts(counts.auto[i])
                counts.cross12[i] = fix_bad_bins_counts(counts.cross12[i])
                counts.cross21[i] = fix_bad_bins_counts(counts.cross21[i])
        kw[name] = counts
    return cls(**kw)


def reshape_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator, n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf, skip_r_bins: int | tuple[int, int] = 0) -> pycorr.twopoint_estimator.BaseTwoPointEstimator:
    """
    Reshape a pycorr estimator to match the settings of the covariance matrix, by rebinning in r and mu, applying r_max cut and skipping some r bins if requested.
    Also checks negative counts and wraps to |mu| bins if the input extends to negative mu.

    Parameters
    ----------
    xi_estimator : pycorr.twopoint_estimator.BaseTwoPointEstimator
        The input pycorr estimator to reshape. It is recommended not to wrap it to |mu| bins beforehand, as the results could be weird if the wrapping was done before rebinning.

    n_mu : integer
        (Optional) the desired number of angular (|mu|) bins (after wrapping to absolute value of mu). If not provided, zero on None, the original number of mu bins is kept.

    r_step : float
        (Optional) the desired width of the radial (separation) bins. If not provided, zero or None, the original radial binning is kept. Rebinning is only supported for linear bins and will be done by integer factor, so the provided r_step needs to be close enough to an integer multiple of the original r_step in ``pycorr``.

    r_max : float
        (Optional) cut the radial bins with separations larger than the given r_max value. If not provided or infinity, no bins are removed by this criterion. The cut is applied after skipping bins but before rebinning.
    
    skip_r_bins : integer or tuple of two integers
        (Optional) removal of some radial bins after radial rebinning.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.
    
    output_cov_file1, output_cov_file2 : string or None
        (Optional) if provided, the text covariance matrices for the corresponding region will be saved in this file.

    Returns
    -------
    xi_estimator_reshaped : pycorr.twopoint_estimator.BaseTwoPointEstimator
        The reshaped pycorr estimator.
    """
    # determine the radius step in pycorr
    if not r_step: r_factor = 1
    else:
        r_steps_orig = np.diff(xi_estimator.edges[0])
        r_step_orig = np.mean(r_steps_orig)
        if not np.allclose(r_steps_orig, r_step_orig, rtol=5e-3, atol=5e-3): raise ValueError("Radial rebinning only supported for linear bins")
        r_factor_exact = r_step / r_step_orig
        r_factor = int(np.rint(r_factor_exact))
        if not np.allclose(r_factor, r_factor_exact, rtol=5e-3): raise ValueError(f"Radial rebinning seems impossible: exact ratio of steps is {r_factor_exact}, closest integer is {r_factor} and that is too far")

    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    # Skip the first bins as requested
    xi_estimator = xi_estimator[skip_r_bins_start * r_factor:]
    # Skip the last bins if requested
    if skip_r_bins_end != 0: xi_estimator = xi_estimator[:-skip_r_bins_end * r_factor] # the slice works unless skip_r_bins_end=0

    # Apply r_max cut
    r_values = xi_estimator.sepavg(axis = 0)
    r_bins_indices = np.where(r_values <= r_max)[0]
    xi_estimator = xi_estimator[:r_bins_indices.max() + 1] # can not apply mask to pycorr estimators; the bins are adjacent anyway

    # determine the mu binning factor
    n_mu_orig = xi_estimator.shape[1]
    need_wrap = np.any(xi_estimator.edges[1] < 0) # wrap if extends to negative mu
    if need_wrap:
        if n_mu_orig % 2 != 0: raise ValueError("Wrapping not possible")
        n_mu_orig //= 2
    if n_mu:
        if n_mu_orig % n_mu != 0: raise ValueError("Angular rebinning not possible")
        mu_factor = n_mu_orig // n_mu
    else: mu_factor = 1 # leave the original number of mu bins

    if need_wrap: xi_estimator = fix_bad_bins_pycorr(xi_estimator) # try to fix bad bins using symmetry
    xi_estimator = xi_estimator[::r_factor, ::mu_factor] # rebin
    if need_wrap: xi_estimator = xi_estimator.wrap() # wrap after rebinning, otherwise the results used to be weird
    elif r_factor != 1 or mu_factor != 1: warn("Seems like the estimator has been wrapped before, nontrivial rebinning after wrapping might be troublesome")

    return xi_estimator
