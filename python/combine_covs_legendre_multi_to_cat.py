"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two full 2-tracer covs following NS/GCcomb procedure for 2 tracers in Legendre mode and convert them . Covariance of N(GC) and S(GC) 2PCF is neglected."

from pycorr import TwoPointCorrelationFunction
import numpy as np
import sys
from .utils import reshape_pycorr
from .convert_cov_legendre_multi import get_cov_header, load_cov_legendre_multi
from .convert_counts_from_pycorr import get_counts_from_pycorr
from .mu_bin_legendre_factors import compute_mu_bin_legendre_factors


def combine_covs_legendre_multi_to_cat(rascalc_results1: str, rascalc_results2: str, pycorr_files1: list[str], pycorr_files2: list[str], output_cov_file: str, max_l: int, r_step: float = 1, skip_r_bins: int = 0, bias1: float = 1, bias2: float = 1, output_cov_file1: str | None = None, output_cov_file2: str | None = None, print_function = print):
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov_legendre_multi(rascalc_results1, max_l, print_function)
    n_bins = len(cov1)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov_legendre_multi(rascalc_results2, max_l, print_function)

    # Read pycorr files to figure out weights
    results1 = []
    for pycorr_file1 in pycorr_files1:
        results1.append(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file1), n_mu_bins = None, r_step = r_step, skip_r_bins = skip_r_bins))
    assert len(results1) == 3, "Wrong number of xi1s"

    n = results1[0].shape[0]
    mu_edges = results1[0].edges[1]

    results2 = []
    for pycorr_file2 in pycorr_files2:
        results2.append(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file2), n_mu_bins = None, r_step = r_step, skip_r_bins = skip_r_bins))
    assert len(results2) == 3, "Wrong number of xi2s"

    # Add weighting by bias for each tracer
    bias_weights = (bias1**2, 2*bias1*bias2, bias2**2) # auto1, cross12, auto2 are multiplied by product of biases of tracers involved in each. Moreover, cross12 enters twice because wrapped cross21 is the same.
    # Now multiply ALL the counts
    results1 = [result * factor for result, factor in zip(results1, bias_weights)]
    results2 = [result * factor for result, factor in zip(results2, bias_weights)]
    # Extract not normalized RR counts
    weights1 = np.array((get_counts_from_pycorr(result, counts_factor = 1) for result in results1))
    weights2 = np.array((get_counts_from_pycorr(result, counts_factor = 1) for result in results2))
    # Other weights will be counts after summation and and normalization
    weight1 = get_counts_from_pycorr(sum(results1).normalize(), counts_factor = 1)
    weight2 = get_counts_from_pycorr(sum(results2).normalize(), counts_factor = 1)

    # Normalize weights
    sum_weight = weight1 + weight2
    weight1 /= sum_weight
    weight2 /= sum_weight
    # Normalize the full weights across correlation function labels
    weights1 /= np.sum(weights1, axis=0)[None, :, :]
    weights2 /= np.sum(weights2, axis=0)[None, :, :]

    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse = True)

    # First, convert multi-tracer cov to single-tracer in each region

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weights1/2[tracer, r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while tracers and radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,tkl,lj,km->tikjm', leg_mu_factors, weights1, mu_leg_factors, np.eye(n)).reshape(n_bins, n_bins // 3)
    pd2 = np.einsum('il,tkl,lj,km->tikjm', leg_mu_factors, weights2, mu_leg_factors, np.eye(n)).reshape(n_bins, n_bins // 3)
    # We have correct [t_in, l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape.
    # The resulting covs are single-tracer (for the combined catalogs) so there is no t_out.

    # Produce single-tracer covs for each region
    cov1 = pd1.T.dot(cov1).dot(pd1)
    cov2 = pd2.T.dot(cov2).dot(pd2)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header = header1) # includes shot-noise rescaling value in the header
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header = header2) # includes shot-noise rescaling value in the header
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    n_bins //= 3 # all covariances are single tracer now

    # Now, combine single-tracer covs

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weight1/2[r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight1, mu_leg_factors, np.eye(n)).reshape(n_bins, n_bins)
    pd2 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight2, mu_leg_factors, np.eye(n)).reshape(n_bins, n_bins)
    # We have correct [l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape

    # Produce and save combined cov
    cov = pd1.T.dot(cov1).dot(pd1) + pd2.T.dot(cov2).dot(pd2)
    np.savetxt(output_cov_file, cov, header=header) # includes source parts and their shot-noise rescaling values in the header

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) not in (13, 15, 17):
        print("Usage: python combine_covs_multi_to_cat.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1_11} {PYCORR_FILE2_11} {PYCORR_FILE1_12} {PYCORR_FILE2_12} {PYCORR_FILE1_22} {PYCORR_FILE2_22} {R_STEP} {MAX_L} {R_BINS_SKIP} {OUTPUT_COV_FILE} [{BIAS1} {BIAS2} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}]].")
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
    bias1 = get_arg_safe(13, float, 1)
    bias2 = get_arg_safe(14, float, 1)
    output_cov_file1 = get_arg_safe(15, str, None)
    output_cov_file2 = get_arg_safe(16, str, None)
    
    combine_covs_legendre_multi_to_cat(rascalc_results1, rascalc_results2, pycorr_files1, pycorr_files2, output_cov_file, max_l, r_step, skip_r_bins, bias1, bias2, output_cov_file1, output_cov_file2)