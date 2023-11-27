"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two full 2-tracer covs following NS/GCcomb procedure for 2 tracers in Legendre mode. Covariance of N(GC) and S(GC) 2PCF is neglected."

from pycorr import TwoPointCorrelationFunction
import numpy as np
import sys
from .utils import reshape_pycorr
from .convert_cov_legendre_multi import get_cov_header, load_cov_legendre_multi
from .convert_counts_from_pycorr import get_counts_from_pycorr
from .mu_bin_legendre_factors import compute_mu_bin_legendre_factors


def combine_covs_legendre_multi(rascalc_results1: str, rascalc_results2: str, pycorr_files1: list[str], pycorr_files2: list[str], output_cov_file: str, max_l: int, r_step: float = 1, skip_r_bins: int = 0, output_cov_file1: str | None = None, output_cov_file2: str | None = None, print_function = print):
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov_legendre_multi(rascalc_results1, max_l, print_function)
    n_bins = len(cov1)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov_legendre_multi(rascalc_results2, max_l, print_function)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header = header1)
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header = header2)
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    # Read pycorr files to figure out weights
    weight1 = []
    for pycorr_file1 in pycorr_files1:
        xi_estimator1 = reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file1), n_mu_bins = None, r_step = r_step, skip_r_bins = skip_r_bins)
        weight1.append(get_counts_from_pycorr(xi_estimator1, counts_factor = 1))
    weight1 = np.array(weight1)

    n_r_bins = xi_estimator1.shape[0]
    mu_edges = xi_estimator1.edges[1]

    weight2 = []
    for pycorr_file2 in pycorr_files2:
        weight2.append(get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file2), n_mu_bins = None, r_step = r_step, skip_r_bins = skip_r_bins), counts_factor = 1))
    weight2 = np.array(weight2)

    # Normalize weights
    sum_weight = weight1 + weight2
    weight1 /= sum_weight
    weight2 /= sum_weight

    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse = True)

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weight1/2[tracer, r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while tracers and radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,tkl,lj,km,tr->tikrjm', leg_mu_factors, weight1, mu_leg_factors, np.eye(n_r_bins), np.eye(3)).reshape(n_bins, n_bins)
    pd2 = np.einsum('il,tkl,lj,km,tr->tikrjm', leg_mu_factors, weight2, mu_leg_factors, np.eye(n_r_bins), np.eye(3)).reshape(n_bins, n_bins)
    # We have correct [t_in, l_in, r_in, t_out, l_out, r_out] ordering and want to make these matrices in the end thus the reshape

    # Produce and save combined cov
    cov = pd1.T.dot(cov1).dot(pd1) + pd2.T.dot(cov2).dot(pd2)
    np.savetxt(output_cov_file, cov, header=header) # includes source parts and their shot-noise rescaling values in the header

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) not in (13, 15):
        print("Usage: python combine_covs_legendre_multi.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1_11} {PYCORR_FILE2_11} {PYCORR_FILE1_12} {PYCORR_FILE2_12} {PYCORR_FILE1_22} {PYCORR_FILE2_22} {N_R_BINS} {MAX_L} {SKIP_R_BINS} {OUTPUT_COV_FILE} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}].")
        sys.exit(1)
    rascalc_results1 = str(sys.argv[1])
    rascalc_results2 = str(sys.argv[2])
    pycorr_files1 = [str(sys.argv[i]) for i in (3, 5, 7)]
    pycorr_files2 = [str(sys.argv[i]) for i in (4, 6, 8)]
    r_step = float(sys.argv[9])
    max_l = int(sys.argv[10])
    skip_r_bins = int(sys.argv[11])
    output_cov_file = str(sys.argv[12])
    from .utils import get_arg_safe
    output_cov_file1 = get_arg_safe(13, str, None)
    output_cov_file2 = get_arg_safe(14, str, None)
    
    combine_covs_legendre_multi(rascalc_results1, rascalc_results2, pycorr_files1, pycorr_files2, output_cov_file, max_l, r_step, skip_r_bins, output_cov_file1, output_cov_file2)