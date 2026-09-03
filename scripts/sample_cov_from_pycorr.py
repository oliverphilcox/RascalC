"This reads cosmodesi/pycorr .npy file(s) and generates sample covariance of xi(s,mu) in text format"

import argparse

parser = argparse.ArgumentParser(description="This script reads cosmodesi/pycorr .npy file(s) and generates sample covariance of xi(s,mu) in text format")
parser.add_argument("pycorr_files", type=str, nargs='+', help="pycorr .npy filenames. in multi-tracer case, you should give all the realizations for the first tracer's auto-correlation, then the same realizations in the same order for the cross-correlation between tracers 1 and 2, and then for the auto-correlation of the second tracer (this can be continued further, but RascalC can not compute semi-analytic covariances for more than 2 tracers at a time.)")
parser.add_argument("output_file", type=str, help="output text file for the sample covariance matrix")
parser.add_argument("r_step", type=float, help="desired width of radial bins")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("r_max", type=float, help="maximum radius (cutoff)")
parser.add_argument("n_correlations", type=int, help="number of correlation functions (1 for single tracer, typically 3 for two tracers)")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pycorr_utils.sample_cov import sample_cov_from_pycorr_files

n_corr: int = args.n_correlations
pycorr_files: list[str] = args.pycorr_files

assert n_corr >= 1, "Need to have at least one correlation"
n_files = len(pycorr_files)
assert n_files % n_corr == 0, "Need to have the same number of files for all correlations"
n_samples = n_files // n_corr
pycorr_files = [pycorr_files[n_samples * i, n_samples * (i+1)] for i in range(n_corr)]

sample_cov_from_pycorr_files(pycorr_files, args.output_file, args.n_mu_bins, args.r_step, args.r_max)