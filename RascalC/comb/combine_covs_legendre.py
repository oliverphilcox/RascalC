"This reads two sets of RascalC results and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure in Legendre mode. Covariance of N and S 2PCF is neglected."

from pycorr import TwoPointCorrelationFunction
import numpy as np
from ..pycorr_utils.utils import reshape_pycorr
from ..cov_utils import get_cov_header, load_cov_legendre
from ..pycorr_utils.counts import get_counts_from_pycorr
from ..mu_bin_legendre_factors import compute_mu_bin_legendre_factors


def combine_covs_legendre(rascalc_results1: str, rascalc_results2: str, pycorr_file1: str, pycorr_file2: str, output_cov_file: str, max_l: int, r_step: float = 1, skip_r_bins: int | tuple[int, int] = 0, output_cov_file1: str | None = None, output_cov_file2: str | None = None, print_function = print):
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov_legendre(rascalc_results1, max_l, print_function)
    n_bins = len(cov1)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov_legendre(rascalc_results2, max_l, print_function)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header = header1)
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header = header2)
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    # Read pycorr files to figure out weights of s, mu binned 2PCF
    xi_estimator1 = reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file1), n_mu = None, r_step = r_step, skip_r_bins = skip_r_bins).normalize()
    n_r_bins = xi_estimator1.shape[0]
    mu_edges = xi_estimator1.edges[1]
    weight1 = get_counts_from_pycorr(xi_estimator1, counts_factor = 1)
    weight2 = get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file2).normalize(), n_mu = None, r_step = r_step, skip_r_bins = skip_r_bins), counts_factor = 1)

    # Normalize weights
    sum_weight = weight1 + weight2
    weight1 /= sum_weight
    weight2 /= sum_weight

    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse = True)

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weight1/2[r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight1, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins)
    pd2 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight2, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins)
    # We have correct [l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape

    # Produce and save combined cov
    cov = pd1.T.dot(cov1).dot(pd1) + pd2.T.dot(cov2).dot(pd2)
    np.savetxt(output_cov_file, cov, header=header) # includes source parts and their shot-noise rescaling values in the header
