## Script to post-process the multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameters, alpha_i, from data derived covariance matrices.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices.
## The effective number of samples, N_eff, is also computed.

import argparse

parser = argparse.ArgumentParser(description="Script to post-process the single-field integrals computed by the C++ code (in radial and angular/mu bins), optimizing the shot-noise rescaling parameter to a data-derived (jackknife) covariance matrix.")
parser.add_argument("xi_jack_files", type=str, nargs=3, help="name of the 3 text files with the jackknife correlation function estimates (in proper RascalC format) in the following order: the first tracer auto-correlation, the cross-correlation of the two tracers, the second tracer auto-correlation")
parser.add_argument("weight_dir", type=str, help="directory containing the jackknife region weights and RR counts")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("output_dir", type=str, help="directory to write the post-processing results")
parser.add_argument("skip_r_bins", type=int, default=0, nargs='?', help="number of last radial/separation bins to discard")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.post_process import post_process_jackknife_multi

post_process_jackknife_multi(args.xi_jack_files[0], args.xi_jack_files[1], args.xi_jack_files[2], args.weight_dir, args.covariance_dir, args.n_mu_bins, args.output_dir, args.skip_r_bins)
