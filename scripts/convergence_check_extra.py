## Script to perform an extra convergence check on full integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Should work in any case - default, jackknife, Legendre, multi-tracer - as it utilizes universal data from RascalC file

import argparse

parser = argparse.ArgumentParser(description="Script to perform an extra convergence check on full integrals; more specifically, divide integral subsamples into halves and check similarity of their average results. Should work in any case - default, jackknife, Legendre, multi-tracer - as it utilizes universal data from a RascalC .npy file.")
parser.add_argument("rascalc_results", type=str, help="RascalC .npy filename")
parser.add_argument("n_subsamples", type=int, default=None, nargs='?', help="number of covariance subsamples to use (optional; by default all are used)")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.convergence_check_extra import convergence_check_extra_file

convergence_check_extra_file(args.rascalc_results, args.n_subsamples, print_function = print)