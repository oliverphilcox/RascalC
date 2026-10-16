## Script to post-process the single-field Legendre binned integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import argparse

parser = argparse.ArgumentParser(description="Script to post-process the single-field Legendre binned integrals computed by the C++ code")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
parser.add_argument("n_r_bins", type=int, help="number of radial/separation bins")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_dir", type=str, help="directory to write the post-processing results")
parser.add_argument("shot_noise_rescaling", type=float, default=1, nargs='?', help="shot-noise rescaling parameter value to use")
parser.add_argument("skip_r_bins", type=int, default=0, nargs='?', help="number of last radial/separation bins to discard")
parser.add_argument("skip_l", type=int, default=0, nargs='?', help="number of last Legendre multipoles to discard")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.post_process import post_process_legendre

post_process_legendre(args.covariance_dir, args.n_r_bins, args.max_l, args.output_dir, args.shot_noise_rescaling, args.skip_r_bins, args.skip_l)