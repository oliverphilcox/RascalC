"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two full 2-tracer covs following NS/GCcomb procedure for 2 tracers in Legendre mode. Covariance of N(GC) and S(GC) 2PCF is neglected."

import argparse

parser = argparse.ArgumentParser(description="This script reads two sets of RascalC results and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure for 2 tracers in Legendre mode. Covariance of N and S 2PCF is neglected.")
parser.add_argument("rascalc_files", type=str, nargs=2, help="RascalC post-processed .npy filenames for the two regions")
parser.add_argument("pycorr_files", type=str, nargs=6, help="6 pycorr .npy filenames for the two regions and 3 distinct correlation functions in the following order: first tracer auto-counts for both regions, then cross-counts between the two tracers for both regions, and finally the second tracer auto-counts for both regions")
parser.add_argument("r_step", type=float, help="desired width of radial bins")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("skip_r_bins", type=int, help="number of last radial/separation bins to discard")
parser.add_argument("output_cov_file", type=str, help="output text file for the covariance in the combined region")
parser.add_argument("output_cov_file1", type=str, default=None, nargs='?', help="optional output text file for the covariance in the first region")
parser.add_argument("output_cov_file2", type=str, default=None, nargs='?', help="optional output text file for the covariance in the second region")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.comb.combine_covs_legendre_multi import combine_covs_legendre_multi

combine_covs_legendre_multi(args.rascalc_files[0], args.rascalc_files[1], [args.pycorr_files[i] for i in range(0, 6, 2)], [args.pycorr_files[i] for i in range(1, 6, 2)], args.output_cov_file, args.max_l, args.r_step, args.skip_r_bins, args.output_cov_file1, args.output_cov_file2)