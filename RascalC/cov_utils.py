"""
This module contains convenience functions that

- read RascalC results (checking eigenvalues of the inversion bias matrix. If they are much smaller than 1, it is safe to simply invert the covariance matrix. Otherwise, a correction factor is necessary.)
- convert (reorder) covariance matrices (as is often needed for Legendre moments)
- and/or export (save) the full covariance matrix to a text file.
"""

import numpy as np
from typing import Callable
from .get_shot_noise_rescaling import get_shot_noise_rescaling


def get_cov_header(rascalc_results_file: str) -> str:
    """
    Format the header for the covariance matrix text file in a consistent way to be used in other functions.
    Currently, just gets the shot-noise rescaling value.
    """
    return "shot_noise_rescaling = " + str(get_shot_noise_rescaling(rascalc_results_file))


def convert_cov_legendre(cov: np.typing.NDArray[np.float64], max_l: int) -> np.typing.NDArray[np.float64]:
    """
    Change the bin ordering of the covariance matrix in Legendre mode for a single tracer.
    The original bin ordering (in ``RascalC`` ``.npy`` files) is by radial bins (top-level) and then by multipoles.
    The resulting bin ordering (for text files) is by multipoles (top-level) and then by radial bins.
    """
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    n_bins = len(cov)
    if n_bins % n_l != 0: raise ValueError("Number of bins in the Legendre covariance must be divisible by the number of even multipoles")
    n_r_bins = n_bins // n_l
    cov = cov.reshape(n_r_bins, n_l, n_r_bins, n_l) # convert to 4D from 2D with [r, l] ordering for both rows and columns
    cov = cov.transpose(1, 0, 3, 2) # change ordering to [l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
    return cov


def convert_cov_legendre_multi(cov: np.typing.NDArray[np.float64], max_l: int) -> np.typing.NDArray[np.float64]:
    """
    Change the bin ordering of the covariance matrix in Legendre mode for two tracers.
    The original bin ordering (in ``RascalC`` ``.npy`` files) is by the correlation function (top-level; first auto-correlation, then cross-correlation and last the second auto-correlation), then by radial bins and then by multipoles (bottom-level).
    The resulting bin ordering (for text files) is by the correlation function (top-level), then by multipoles and then by radial bins (bottom-level).
    """
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    n_bins = len(cov)
    if n_bins % (3 * n_l) != 0: raise ValueError("Number of bins in the multi-tracer Legendre covariance must be divisible by thrice the number of even multipoles")
    n_r_bins = n_bins // (3 * n_l)
    cov = cov.reshape(3, n_r_bins, n_l, 3, n_r_bins, n_l) # convert to 6D from 2D with [t, r, l] ordering for both rows and columns
    cov = cov.transpose(0, 2, 1, 3, 5, 4) # change ordering to [t, l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 6D to 2D
    return cov


def convert_cov_multi_to_cross(cov: np.typing.NDArray[np.float64]) -> np.typing.NDArray[np.float64]:
    """
    Select only the cross x cross-correlation part (middle block) of the full two-tracer covariance.
    This function does not change the bin order of the covariance.
    """
    n_bins = len(cov)
    if n_bins % 3 != 0: raise ValueError("Number of bins in the multi-tracer covariance must be divisible by 3")
    n_bins //= 3
    cov = cov.reshape(3, n_bins, 3, n_bins) # convert to 4D from 2D with [t, b] ordering for both rows and columns
    cov = cov.transpose(0, 2, 1, 3) # change ordering to [t1, t2, b1, b2]
    cov = cov[1, 1] # select the covariance for cross-correlation only and reduce from 4D to 2D
    return cov


def convert_cov_legendre_multi_to_cross(cov: np.typing.NDArray[np.float64], max_l: int) -> np.typing.NDArray[np.float64]:
    """
    Select only the cross x cross-correlation part (middle block) of the full two-tracer covariance in Legendre mode.
    This function also changes the bin order of the covariance as in :func:`convert_cov_legendre`.
    """
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    n_bins = len(cov)
    if n_bins % (3 * n_l) != 0: raise ValueError("Number of bins in the multi-tracer Legendre covariance must be divisible by thrice the number of even multipoles")
    n_bins //= 3
    n_r_bins = n_bins // n_l
    cov = cov.reshape(3, n_r_bins, n_l, 3, n_r_bins, n_l) # convert to 6D from 2D with [t, r, l] ordering for both rows and columns
    cov = cov.transpose(0, 3, 2, 1, 5, 4) # change ordering to [t1, t2, l1, r1, l2, r2]
    cov = cov[1, 1] # select the covariance for cross-correlation only and reduce from 6D to 4D with [l, r] ordering for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
    return cov


def load_cov(rascalc_results_file: str, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    "Load the theoretical covariance matrix from RascalC results file as-is, intended for the s_mu mode."
    with np.load(rascalc_results_file) as f:
        print_function(f"Max abs eigenvalue of inversion bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
        # if the printed value is small the cov matrix should be safe to invert as is
        return f['full_theory_covariance']


def load_cov_legendre(rascalc_results_file: str, max_l: int, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    "Load the theoretical covariance matrix from RascalC results file and change the bin ordering as in :func:`convert_cov_legendre`; intended for Legendre single-tracer mode."
    return convert_cov_legendre(load_cov(rascalc_results_file, print_function), max_l)


def load_cov_legendre_multi(rascalc_results_file: str, max_l: int, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    "Load the theoretical covariance matrix from RascalC results file and change the bin ordering as in :func:`convert_cov_legendre_multi`; intended for Legendre two-tracer mode."
    return convert_cov_legendre_multi(load_cov(rascalc_results_file, print_function), max_l)


def load_cov_cross_from_multi(rascalc_results_file: str, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    "Load only the cross x cross-correlation part (middle block) of the full two-tracer covariance in a RascalC results file without conversion."
    return convert_cov_multi_to_cross(load_cov(rascalc_results_file, print_function))


def load_cov_legendre_cross_from_multi(rascalc_results_file: str, max_l: int, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    "Load only the cross x cross-correlation part (middle block) of the full two-tracer covariance in a RascalC results file with the Legendre-mode conversion (see :func:`convert_cov_legendre`)."
    return convert_cov_legendre_multi_to_cross(load_cov(rascalc_results_file, print_function), max_l)


def export_cov(rascalc_results_file: str, output_cov_file: str, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix from RascalC results file to a text file as-is, intended for the s_mu mode."
    np.savetxt(output_cov_file, load_cov(rascalc_results_file, print_function = print_function), header = get_cov_header(rascalc_results_file))


def export_cov_legendre(rascalc_results_file: str, max_l: int, output_cov_file: str, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix from RascalC results file to a text file with conversion appropriate for single-tracer Legendre modes (see :func:`convert_cov_legendre`)."
    np.savetxt(output_cov_file, load_cov_legendre(rascalc_results_file, max_l, print_function = print_function), header = get_cov_header(rascalc_results_file))


def export_cov_legendre_multi(rascalc_results_file: str, max_l: int, output_cov_file: str, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix from RascalC results file to a text file with conversion appropriate for two-tracer Legendre modes (see :func:`convert_cov_legendre_multi`)."
    np.savetxt(output_cov_file, load_cov_legendre_multi(rascalc_results_file, max_l, print_function = print_function), header = get_cov_header(rascalc_results_file))


def export_cov_cross(rascalc_results_file: str, output_cov_file: str, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix for cross-correlation only from RascalC results file to a text file without conversion, intended for the s_mu mode."
    np.savetxt(output_cov_file, load_cov_cross_from_multi(rascalc_results_file, print_function = print_function), header = get_cov_header(rascalc_results_file))


def export_cov_legendre_cross(rascalc_results_file: str, max_l: int, output_cov_file: str, print_function: Callable[[str], None] = print) -> None:
    "Export the theoretical covariance matrix for cross-correlation only from RascalC results file to a text file with conversion appropriate for the Legendre modes (see :func:`convert_cov_legendre`)."
    np.savetxt(output_cov_file, load_cov_legendre_cross_from_multi(rascalc_results_file, max_l, print_function = print_function), header = get_cov_header(rascalc_results_file))


def convert_txt_cov_multi_to_cross(input_file: str, output_file: str) -> None:
    "Convert a plain-text covariance file from full two-tracer to cross x cross-correlation only by selecting the middle block of the matrix. Should be appropriate for any mode."
    with open(input_file) as f:
        header = f.readline()
    if header.startswith("#"): header = header.strip("#").strip() # first remove the comment character and then any whitespaces
    else: header = None
    np.savetxt(output_file, convert_cov_multi_to_cross(np.loadtxt(input_file)), header = "converted from " + input_file + (", " + header) * bool(header))
