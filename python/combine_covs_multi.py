"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure for 2 tracers. Covariance of N and S 2PCF is neglected."

from pycorr import TwoPointCorrelationFunction
import numpy as np
import sys
from utils import reshape_pycorr
from convert_cov import get_cov_header, load_cov
from convert_counts_from_pycorr import get_counts_from_pycorr

def combine_covs_multi(rascalc_results1: str, rascalc_results2: str, pycorr_files1: list[str], pycorr_files2: list[str], output_cov_file: str, n_mu_bins: int | None = None, r_step: float = 1, skip_r_bins: int = 0, output_cov_file1: str | None = None, output_cov_file2: str | None = None, print_function = print):
    if len(pycorr_files1) != len(pycorr_files2): raise ValueError("Need the same number of pycorr files for both results")
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
    weight1 = np.zeros(0)
    for pycorr_file1 in pycorr_files1:
        weight1 = np.append(weight1, get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file1), n_mu_bins, r_step, skip_r_bins = skip_r_bins), counts_factor = 1).ravel())
    weight2 = np.zeros(0)
    for pycorr_file2 in pycorr_files2:
        weight2 = np.append(weight2, get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file2), n_mu_bins, r_step, skip_r_bins = skip_r_bins), counts_factor = 1).ravel())

    # Produce and save combined cov
    # following xi = (xi1 * weight1 + xi2 * weight2) / (weight1 + weight2)
    cov = (cov1 * weight1[None, :] * weight1[:, None] + cov2 * weight2[None, :] * weight2[:, None]) / (weight1 + weight2)[None, :] / (weight1 + weight2)[:, None]
    np.savetxt(output_cov_file, cov, header = header) # includes source parts and their shot-noise rescaling values in the header

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) not in (13, 15):
        print("Usage: python combine_covs_multi.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1_11} {PYCORR_FILE2_11} {PYCORR_FILE1_12} {PYCORR_FILE2_12} {PYCORR_FILE1_22} {PYCORR_FILE2_22} {R_STEP} {N_MU_BINS} {SKIP_R_BINS} {OUTPUT_COV_FILE} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}].")
        sys.exit(1)
    rascalc_results1 = str(sys.argv[1])
    rascalc_results2 = str(sys.argv[2])
    pycorr_files1 = [str(sys.argv[i]) for i in (3, 5, 7)]
    pycorr_files2 = [str(sys.argv[i]) for i in (4, 6, 8)]
    r_step = float(sys.argv[9])
    n_mu_bins = int(sys.argv[10])
    skip_r_bins = int(sys.argv[11])
    output_cov_file = str(sys.argv[12])
    from utils import get_arg_safe
    output_cov_file1 = get_arg_safe(13, str, None)
    output_cov_file2 = get_arg_safe(14, str, None)
    
    combine_covs_multi(rascalc_results1, rascalc_results2, pycorr_files1, pycorr_files2, output_cov_file, n_mu_bins, r_step, skip_r_bins, output_cov_file1, output_cov_file2)