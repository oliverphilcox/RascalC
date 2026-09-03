"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC Legendre mode results and saves full cov to text file, changing the indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix"

import argparse

parser = argparse.ArgumentParser(description="Simple convenience script that reads RascalC results and saves full cov to text file, changing the single-tracer Legendre-mode indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix")
parser.add_argument("rascalc_file", type=str, help="RascalC post-processed .npy filename")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_cov_file", type=str, help="output text file for the covariance")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.cov_utils import export_cov_legendre

export_cov_legendre(args.rascalc_file, args.max_l, args.output_cov_file)
