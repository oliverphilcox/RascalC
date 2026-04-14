import lsstypes
import numpy as np
from warnings import warn


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
    # fixes bins with negative wcounts by overwriting their content by reflection
    # only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32
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
