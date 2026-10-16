## Script to post-process the multi-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import argparse

parser = argparse.ArgumentParser(description="Script to post-process the multi-field integrals computed by the C++ code (in radial and angular/mu bins).")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
parser.add_argument("n_r_bins", type=int, help="number of radial/separation bins")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("output_dir", type=str, help="directory to write the post-processing results")
parser.add_argument("shot_noise_rescaling1", type=float, default=1, nargs='?', help="shot-noise rescaling parameter value to use for the first tracer")
parser.add_argument("shot_noise_rescaling2", type=float, default=1, nargs='?', help="shot-noise rescaling parameter value to use for the second tracer")
parser.add_argument("skip_r_bins", type=int, default=0, nargs='?', help="number of last radial/separation bins to discard")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.post_process import post_process_default_multi

post_process_default_multi(args.covariance_dir, args.n_r_bins, args.n_mu_bins, args.output_dir, args.shot_noise_rescaling1, args.shot_noise_rescaling2, args.skip_r_bins)
