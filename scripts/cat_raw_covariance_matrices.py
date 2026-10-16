## Script to catenate (take a subset of subsamples or resample) raw covariance matrix files
## More specifically, copy part of partial results to other directory and recompute totals by averaging
## Determines single-field vs multi-field and jackknife automatically
## Do not use if subsamples have different numbers of pairs/triples/quadruplets

import argparse

parser = argparse.ArgumentParser(description="Script to catenate (take a subset of subsamples or resample) raw covariance matrix files or directories. More specifically, copy part of partial results to other directory and recompute totals by averaging. Determines single-field vs multi-field and jackknife automatically. Do not use if subsamples have different numbers of pairs/triples/quadruplets.")
parser.add_argument("n_r_bins", type=int, help="number of radial/separation bins")
parser.add_argument("n_mu_bins_or_max_l_specifier", type=int, help="string specifying the number of angular/mu bins or max_l (Legendre), the format is m{n_mu_bins} or l{max_l} accordingly")
parser.add_argument("input_covariance_dirs", type=str, nargs='+', help="one directory or multiple directories to read the covariance matrix subsamples from")
parser.add_argument("--ns_subsamples", type=int, nargs='*', default=None, help="number of subsamples to use, one number for each covariance directory")
parser.add_argument("--collapse_factor", type=int, default=1, help="reduce the number of samples by this factor by averaging batches of this size (default 1 = no reduction)")
parser.add_argument("output_covariance_dir", type=str, help="directory to write the resulting covariance matrices and subsamples")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.raw_covariance_matrices import cat_raw_covariance_matrices

cat_raw_covariance_matrices(args.n, args.mstr, args.input_covariance_dirs, args.ns_subsamples or [None] * len(args.input_covariance_dirs), args.output_covariance_dir, args.collapse_factor)