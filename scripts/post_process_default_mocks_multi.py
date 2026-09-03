## Script to post-process the multi-field integrals computed by the C++ code. This computes two shot-noise rescaling parameters, alphas, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import argparse

parser = argparse.ArgumentParser(description="Script to post-process the multi-field integrals computed by the C++ code (in radial and angular/mu bins), optimizing the two shot-noise rescaling parameters to the mock-derived covariance matrix.")
parser.add_argument("mock_cov_file", type=str, help="name of the text file with the mock sample covariance matrix")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
parser.add_argument("n_r_bins", type=int, help="number of radial/separation bins")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("output_dir", type=str, help="directory to write the post-processing results")
parser.add_argument("skip_r_bins", type=int, default=0, nargs='?', help="number of last radial/separation bins to discard")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.post_process import post_process_default_mocks_multi

post_process_default_mocks_multi(args.mock_cov_file, args.covariance_dir, args.n_r_bins, args.n_mu_bins, args.output_dir, args.skip_r_bins)
