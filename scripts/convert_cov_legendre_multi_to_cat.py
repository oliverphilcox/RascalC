"This reads a .npy file of RascalC Legendre results (or txt covariance converted previously) and a triplet of cosmodesi/pycorr .npy files to produce a covariance for a catalog of these two tracers concatenated."

import argparse

parser = argparse.ArgumentParser(description="This reads a .npy file of RascalC Legendre results (or txt covariance converted previously) and a triplet of cosmodesi/pycorr .npy files to produce a covariance for a catalog of these two tracers concatenated (each tracer upweighted by its linear bias).")
parser.add_argument("rascalc_file", type=str, help="RascalC post-processed .npy filename")
parser.add_argument("pycorr_files", type=str, nargs=3, help="3 pycorr .npy filenames for the distinct correlation functions in the following order: first tracer auto-counts, then cross-counts between the two tracers, and finally the second tracer auto-counts")
parser.add_argument("r_step", type=float, help="desired width of radial bins")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("skip_r_bins", type=int, help="number of last radial/separation bins to discard")
parser.add_argument("output_cov_file", type=str, help="output text file for the combined tracer covariance")
parser.add_argument("bias1", type=float, default=1, nargs='?', help="linear bias for the first tracer (default 1)")
parser.add_argument("bias2", type=float, default=1, nargs='?', help="linear bias for the second tracer (default 1)")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.comb.convert_cov_legendre_multi_to_cat import convert_cov_legendre_multi_to_cat

convert_cov_legendre_multi_to_cat(args.rascalc_file, args.pycorr_files, args.output_cov_file, args.max_l, args.r_step, args.skip_r_bins, args.bias1, args.bias2)