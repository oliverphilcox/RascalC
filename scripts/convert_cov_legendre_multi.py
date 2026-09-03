"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC 2-tracer Legendre mode results and saves full cov to text file, changing the indexing from [t, r, l] to [t, l, r], in addition checking eigenvalues of bias matrix"

import argparse

parser = argparse.ArgumentParser(description="Simple convenience script that reads RascalC results and saves full cov to text file, changing the 2-tracer Legendre-mode indexing from [t, r, l] to [t, l, r], in addition checking eigenvalues of bias matrix")
parser.add_argument("rascalc_file", type=str, help="RascalC post-processed .npy filename")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_cov_file", type=str, help="output text file for the covariance")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.cov_utils import export_cov_legendre_multi

export_cov_legendre_multi(args.rascalc_file, args.max_l, args.output_cov_file)
