"""
This module contains convenience functions that

- read RascalC results (checking eigenvalues of the inversion bias matrix. If they are much smaller than 1, it is safe to simply invert the covariance matrix. Otherwise, a correction factor is necessary.
- convert (reorder) 3PCF covariance matrices (as is usually needed for Legendre moments)
- and/or export (save) the full 3PCF covariance matrix to a text file.
"""

import numpy as np
import numpy.typing as npt
from typing import Callable
from .get_shot_noise_rescaling import get_shot_noise_rescaling # for convenience
from .cov_utils import get_cov_header, load_cov


def convert_cov_3pcf_legendre(cov: npt.NDArray[np.float64], max_l: int, exclude_odd_l: bool = False, apply_scaling: bool = True) -> npt.NDArray[np.float64]:
    """
    Change the bin ordering of the 3PCF covariance matrix in Legendre mode for a single tracer.
    The bin ordering in ``RascalC`` ``.npy`` files is by radial bin pairs (top-level) and then by multipoles.
    The resulting bin ordering (for text files) is by multipoles (top-level) and then by radial bin pairs.
    By default, the ell-dependent scaling factor is applied to convert to ENCORE basis functions from RascalC Legendre polynomial basis, this can be disabled with ``apply_scaling=False``.
    """
    n_l = max_l // (1 + exclude_odd_l) + 1
    n_bins = len(cov)
    if n_bins % n_l != 0: raise ValueError("Number of bins in the Legendre covariance must be divisible by the number of multipoles (even and odd)")
    n_r_bin_pairs = n_bins // n_l
    cov = cov.reshape(n_r_bin_pairs, n_l, n_r_bin_pairs, n_l) # convert to 4D from 2D with [r_bin_pair, l] ordering for both rows and columns
    if apply_scaling: # prepare to apply the ell-dependent factor between RascalC and ENCORE basis functions
        ells = np.arange(0, max_l + 1, 1 + exclude_odd_l)
        ell_factor = ((-1)**ells * np.sqrt(2 * ells + 1) / (4 * np.pi)) # the ell-dependent factor between the ENCORE 3-point basis functions and Legendre polynomials given by Equation (17) in https://arxiv.org/pdf/2105.08722
        cov /= ell_factor[:, None, None] * ell_factor[None, None, :] # apply the factor, NumPy broadcasting matches trailing dimensions. need to check if it is not multiplication; there might also be a factor of 2 or something similar
    cov = cov.transpose(1, 0, 3, 2) # change ordering to [l, r_pair] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
    return cov


def load_cov_3pcf_legendre(rascalc_results_file: str, max_l: int, exclude_odd_l: bool = False, apply_scaling: bool = True, print_function: Callable[[str], None] = print) -> npt.NDArray[np.float64]:
    "Load the theoretical covariance matrix from RascalC results file, change the bin ordering and apply scaling as in :func:`convert_cov_3pcf_legendre`; intended for 3PCF Legendre single-tracer mode."
    return convert_cov_3pcf_legendre(load_cov(rascalc_results_file, print_function), max_l, exclude_odd_l=exclude_odd_l, apply_scaling=apply_scaling)


def export_cov_3pcf_legendre(rascalc_results_file: str, max_l: int, output_cov_file: str, exclude_odd_l: bool = False, apply_scaling: bool = True, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix from RascalC results file to a text file with conversion appropriate for 3PCF single-tracer Legendre mode (see :func:`convert_cov_3pcf_legendre`)."
    np.savetxt(output_cov_file, load_cov_3pcf_legendre(rascalc_results_file, max_l, exclude_odd_l=exclude_odd_l, apply_scaling=apply_scaling, print_function = print_function), header = get_cov_header(rascalc_results_file))