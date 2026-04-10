import lsstypes
import numpy as np


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