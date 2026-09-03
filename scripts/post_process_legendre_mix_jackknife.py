## Script to post-process the single-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import argparse

parser = argparse.ArgumentParser(description="Script to post-process the single-field Legendre binned integrals computed by the C++ code, optimizing the shot-noise rescaling parameter to a data-derived (jackknife) covariance matrix.")
parser.add_argument("xi_jack_file", type=str, help="name of the text file with the jackknife correlation function estimates (in proper RascalC format)")
parser.add_argument("weight_dir", type=str, help="directory containing the jackknife region weights and RR counts")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_dir", type=str, help="directory to write the post-processing results")
parser.add_argument("skip_r_bins", type=int, default=0, nargs='?', help="number of last radial/separation bins to discard")
parser.add_argument("skip_l", type=int, default=0, nargs='?', help="number of last Legendre multipoles to discard")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.post_process import post_process_legendre_mix_jackknife

post_process_legendre_mix_jackknife(args.xi_jack_file, args.weight_dir, args.file_root, args.n_mu_bins, args.max_l, args.output_dir, args.skip_r_bins, args.skip_l)