"This reads two sets of RascalC results and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure. Covariance of N and S 2PCF is neglected."

from pycorr import TwoPointCorrelationFunction
import numpy as np
from ..pycorr_utils.utils import reshape_pycorr
from ..cov_utils import get_cov_header, load_cov
from ..pycorr_utils.counts import get_counts_from_pycorr


def combine_covs(rascalc_results1: str, rascalc_results2: str, pycorr_file1: str, pycorr_file2: str, output_cov_file: str, n_mu_bins: int | None = None, r_step: float = 1, skip_r_bins: int = 0, output_cov_file1: str | None = None, output_cov_file2: str | None = None, print_function = print):
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov(rascalc_results1, print_function)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov(rascalc_results2, print_function)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header = header1)
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header = header2)
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    # Read pycorr files to figure out weights
    weight1 = get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file1), n_mu_bins, r_step, skip_r_bins = skip_r_bins).normalize(), counts_factor = 1).ravel()
    weight2 = get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file2), n_mu_bins, r_step, skip_r_bins = skip_r_bins).normalize(), counts_factor = 1).ravel()

    # Produce and save combined cov
    # following xi = (xi1 * weight1 + xi2 * weight2) / (weight1 + weight2)
    cov = (cov1 * weight1[None, :] * weight1[:, None] + cov2 * weight2[None, :] * weight2[:, None]) / (weight1 + weight2)[None, :] / (weight1 + weight2)[:, None]
    np.savetxt(output_cov_file, cov, header = header) # includes source parts and their shot-noise rescaling values in the header
