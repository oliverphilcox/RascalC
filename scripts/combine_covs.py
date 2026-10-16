"This reads two sets of RascalC results and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure. Covariance of N and S 2PCF is neglected."

import argparse

parser = argparse.ArgumentParser(description="This script reads two sets of RascalC results (in radial/separation and angular/mu bins) and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure. Covariance of N and S 2PCF is neglected.")
parser.add_argument("rascalc_files", type=str, nargs=2, help="RascalC post-processed .npy filenames for the two regions")
parser.add_argument("pycorr_files", type=str, nargs=2, help="pycorr .npy filenames for the two regions")
parser.add_argument("r_step", type=float, help="desired width of radial bins")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins (0 means to preserve the original number)")
parser.add_argument("skip_r_bins", type=int, help="number of last radial/separation bins to discard")
parser.add_argument("output_cov_file", type=str, help="output text file for the covariance in the combined region")
parser.add_argument("output_cov_file1", type=str, default=None, nargs='?', help="optional output text file for the covariance in the first region")
parser.add_argument("output_cov_file2", type=str, default=None, nargs='?', help="optional output text file for the covariance in the second region")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.comb.combine_covs import combine_covs

combine_covs(args.rascalc_files[0], args.rascalc_files[1], args.pycorr_files[0], args.pycorr_files[1], args.output_cov_file, args.n_mu_bins or None, args.r_step, args.skip_r_bins, args.output_cov_file1, args.output_cov_file2)