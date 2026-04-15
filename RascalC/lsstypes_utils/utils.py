import lsstypes
import numpy as np
from warnings import warn
from ..utils import format_skip_r_bins


# copied from https://github.com/adematti/lsstypes/blob/3bf32b393f81fa7068fbccd027fa793193e056c3/tests/test.py#L1044
def wrap_counts(count: lsstypes.Count2) -> lsstypes.Count2:
    """
    Wrap a Count2 with symmetric edges along the 2nd dimension as a Count2 with half the number of bins along that dimension, and counts summed accordingly.
    This is useful for example to wrap a correlation function in (s, mu) coordinates as a correlation function in (s, |mu|) coordinates.
    """
    if count.shape[1] % 2:
        raise ValueError(f'input counts cannot be wrapped as 2nd dimension is odd = {count.shape[1]:d}')
    mid = count.shape[1] // 2
    sl_neg, sl_pos = slice(mid - 1, None, -1), slice(mid, None, 1)
    coord_name = count._coords_names[1]
    mu_edges = count.edges(coord_name)
    if not np.allclose(mu_edges[sl_neg], - mu_edges[sl_pos, ::-1]):
        raise ValueError(f'input counts cannot be wrapped as 2nd dimension edges are not symmetric; {mu_edges[sl_neg]} != {-mu_edges[sl_pos, ::-1]}')
    mu_edges = mu_edges[sl_pos]
    counts = count.values('counts')
    # Sum counts in the negative and positive halves along the 2nd dimension
    counts = counts[..., sl_neg] + counts[..., sl_pos]
    norm = count.values('norm')[..., sl_pos]
    # Prepare wrapped edges
    edges = dict(count.edges())
    edges[coord_name] = mu_edges
    edges = {'{}_edges'.format(coord): edge for coord, edge in edges.items()}
    # Prepare wrapped coordinates
    mu_coord = count.coords(coord_name)[sl_pos]
    coords = dict(count.coords())
    coords[coord_name] = mu_coord
    return lsstypes.Count2(counts=counts, norm=norm, **coords, **edges, coords=list(coords), attrs=dict(count.attrs))


def wrap_correlation(corr: lsstypes.Count2Correlation) -> lsstypes.Count2Correlation:
    """
    Wrap a Count2Correlation with symmetric edges along the 2nd dimension as a Count2Correlation with half the number of bins along that dimension, and counts summed accordingly.
    This is useful for example to wrap a correlation function in (s, mu) coordinates as a correlation function in (s, |mu|) coordinates.
    """
    # copied from https://github.com/adematti/lsstypes/blob/3bf32b393f81fa7068fbccd027fa793193e056c3/tests/test.py#L1077
    return corr.map(lambda count: wrap_counts(count), level=None, is_leaf=lambda branch: False)  #type(branch) is lsstypes.Count2)


def fix_bad_bins_lsstypes(xi_estimator: lsstypes.Count2Correlation) -> lsstypes.Count2Correlation:
    """
    Fixes bins with negative counts by overwriting their content by reflection.
    Only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32.
    """
    kw = {}
    for count_name in xi_estimator.count_names:
        counts : lsstypes.Count2 = xi_estimator.get(count_name)
        bad_bins_mask = counts.values('normalized_counts') < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {count_name}.normalized_counts ({counts.values('normalized_counts')[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.values('normalized_counts')[s_bin, -1-mu_bin]:.2e})")
            counts._data['normalized_counts'][s_bin, mu_bin] = counts._data['normalized_counts'][s_bin, -1-mu_bin]
            counts._data['norm'][s_bin, mu_bin] = counts._data['norm'][s_bin, -1-mu_bin]
        bad_bins_mask = counts.values('counts') < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {count_name}.counts ({counts.values('counts')[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.values('counts')[s_bin, -1-mu_bin]:.2e})")
            counts._data['normalized_counts'][s_bin, mu_bin] = counts._data['normalized_counts'][s_bin, -1-mu_bin]
            counts._data['norm'][s_bin, mu_bin] = counts._data['norm'][s_bin, -1-mu_bin]
        kw[count_name] = counts
    return lsstypes.Count2Correlation(**kw, estimator=xi_estimator.estimator, attrs=xi_estimator.attrs)


def reshape_lsstypes(xi_estimator: lsstypes.Count2Correlation, n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf, skip_r_bins: int | tuple[int, int] = 0) -> lsstypes.Count2Correlation:
    """
    Reshape a lsstypes estimator to match the settings of the covariance matrix, by rebinning in r and mu, applying r_max cut and skipping some r bins if requested.
    Also checks negative counts and wraps to |mu| bins if the input extends to negative mu.

    Parameters
    ----------
    xi_estimator : lsstypes.Count2Correlation
        The input pycorr estimator to reshape. It is recommended not to wrap it to |mu| bins beforehand, as the results could be weird (in pycorr) if the wrapping was done before rebinning.

    n_mu : integer
        (Optional) the desired number of angular (|mu|) bins (after wrapping to absolute value of mu). If not provided, zero on None, the original number of mu bins is kept.

    r_step : float
        (Optional) the desired width of the radial (separation) bins. If not provided, zero or None, the original radial binning is kept. Rebinning is only supported for linear bins and will be done by integer factor, so the provided r_step needs to be close enough to an integer multiple of the original r_step in ``lsstypes``.

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
    xi_estimator_reshaped : lsstypes.Count2Correlation
        The reshaped lsstypes estimator.
    """
    # determine the radius step in lsstypes
    if not r_step: r_factor = 1
    else:
        r_steps_orig = np.diff(xi_estimator.get('RR')._data['s_edges']) # s_edges is a 2D array with shape (n_s_bins, 2), where the second dimension corresponds to the left and right edges of the bins, so diff along axis 0 gives the widths of the bins; we check that they are all close enough to their mean to make sure the bins are linear
        r_step_orig = np.mean(r_steps_orig)
        if not np.allclose(r_steps_orig, r_step_orig, rtol=5e-3, atol=5e-3): raise ValueError("Radial rebinning only supported for linear bins")
        r_factor_exact = r_step / r_step_orig
        r_factor = int(np.rint(r_factor_exact))
        if not np.allclose(r_factor, r_factor_exact, rtol=5e-3): raise ValueError(f"Radial rebinning seems impossible: exact ratio of steps is {r_factor_exact}, closest integer is {r_factor} and that is too far")

    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    # skip the bins as requested
    xi_estimator = xi_estimator.select(s=slice(skip_r_bins_start * r_factor, None if skip_r_bins_end == 0 else -skip_r_bins_end * r_factor))

    # Apply r_max cut
    xi_estimator = xi_estimator.select(s=(0, r_max))

    # determine the mu binning factor
    mu_edges = xi_estimator.get('RR')._data['mu_edges'] # mu_edges is a 2D array with shape (n_mu_bins, 2), where the second dimension corresponds to the left and right edges of the bins
    n_mu_orig = len(mu_edges)
    need_wrap = np.any(mu_edges < 0) # wrap if extends to negative mu
    if need_wrap:
        if n_mu_orig % 2 != 0: raise ValueError(f"Wrapping not possible")
        n_mu_orig //= 2
    if n_mu:
        if n_mu_orig % n_mu != 0: raise ValueError("Angular rebinning not possible")
        mu_factor = n_mu_orig // n_mu
    else: mu_factor = 1 # leave the original number of mu bins

    if need_wrap: xi_estimator = fix_bad_bins_lsstypes(xi_estimator) # try to fix bad bins using symmetry
    xi_estimator = xi_estimator.select(s=slice(0, None, r_factor), mu=slice(0, None, mu_factor)) # rebin
    if need_wrap: xi_estimator = wrap_correlation(xi_estimator) # wrap after rebinning, otherwise the results used to be weird (in pycorr, but let's just leave it so here as well)
    elif r_factor != 1 or mu_factor != 1: warn("Seems like the estimator has been wrapped before, nontrivial rebinning after wrapping might be troublesome (in pycorr)")

    return xi_estimator
