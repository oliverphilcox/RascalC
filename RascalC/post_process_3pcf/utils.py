import numpy as np
from ..utils import symmetrized, format_skip_r_bins
from ..post_process.utils import apply_cov_filter


def cov_filter_3pcf_legendre(n: int, max_l: int, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, exclude_samebins: bool = True, exclude_odd_l: bool = False):
    """Produce a 2D indexing array for 3PCF Legendre covariance matrices in RascalC convention (not yet fully compatible with ENCORE bin ordering)."""
    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    n_l = max_l + 1
    l_indices = np.arange(0, n_l, 1+exclude_odd_l)
    if skip_l > 0: l_indices = l_indices[:-skip_l] # without the condition, wouldn't work right for skip_l=0
    r_indices = np.arange(skip_r_bins_start, n - skip_r_bins_end)
    r_indices1, r_indices2 = [a.ravel() for a in np.meshgrid(r_indices, r_indices, indexing='ij')] # flattened array indices
    r_filter = (r_indices1 <= r_indices2 - exclude_samebins) # strictly less for exclude_samebin=True, less or equal otherwise
    r_indices1, r_indices2 = r_indices1[r_filter], r_indices2[r_filter] # apply filter to both index arrays
    indices_l_r = (n_l * (n * r_indices1 + r_indices2))[:, None] + l_indices[None, :]
    # indices_l_r = (n_l * (n * r_indices1 + r_indices2))[None, :] + l_indices[:, None] # could switch the l vs bin pair ordering right here easily but decided not to
    indices_1d = indices_l_r.ravel()
    return np.ix_(indices_1d, indices_1d)


def symmetrized_3pcf(A: np.typing.NDArray[np.float64], n: int, max_l: int) -> np.typing.NDArray[np.float64]:
    "Symmetrize a 3PCF covariance matrix (2D array), or an array of 3PCF covariance matrices (3D array)"
    if len(A.shape) not in (2, 3): raise ValueError("Dimension of the input array must be 2 or 3")
    m = max_l + 1
    if not np.array_equal(A.shape[-2:], [n * n * m] * 2): raise ValueError("Unexpected shape in the last 2 dimensions")
    leading_dims = list(A.shape[:-2]) # list containing a leading dimension for the array of covariance matrices, and empty for a single covariance matrix
    A1 = A.reshape(leading_dims + [n, n, m] * 2) # last 6 axes will be [r1, r2, l12, r3, r4, l34]
    A2 = (A1 + A1.swapaxes(-2, -3)) / 2 # symmetrize wrt swaps of r3 and r4. Create new array against the risk of A1 being a view of original A
    A2 = (A2 + A2.swapaxes(-5, -6)) / 2 # symmetrize wrt swaps of r1 and r2
    A2 = A2.reshape(leading_dims + [n * n * m] * 2) # back to the original shape
    return symmetrized(A2) # finally, symetrize wrt full covariance matrix bin swap


def load_matrices(input_data: dict[str], n: int, max_l: int, cov_filter: np.typing.NDArray[np.int_], full: bool = True) -> tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64], np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]:
    """Load the 3PCF single-tracer covariance matrix terms."""
    matrices = []
    for npoints in range(3, 7):
        these_matrices = [input_data[f"c{npoints}_{index}" + "_full" * full] for index in range(npoints == 6, 2)] # exclude c6_0 term, because it should be small but is also hard to compute (see Section 5.2.3 and Appendix A of https://arxiv.org/abs/1910.04764)
        this_matrix = these_matrices[0]
        if npoints != 6: this_matrix += these_matrices[1]
        matrices.append(apply_cov_filter(symmetrized_3pcf(this_matrix, n, max_l), cov_filter)) # symmetrize before filtering, because filtering removes repeating bin pairs
    return tuple(matrices)


def add_cov_terms(c3: np.typing.NDArray[np.float64], c4: np.typing.NDArray[np.float64], c5: np.typing.NDArray[np.float64], c6: np.typing.NDArray[np.float64], alpha: float = 1) -> np.typing.NDArray[np.float64]:
    """Add the 3PCF single-tracer covariance matrix terms with a given shot-noise rescaling value."""
    return c6 + c5 * alpha + c4 * alpha**2 + c3 * alpha**3