"Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of the bias matrix"

import argparse

parser = argparse.ArgumentParser(description="Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of the bias matrix")
parser.add_argument("rascalc_file", type=str, help="RascalC post-processed .npy filename")
parser.add_argument("output_cov_file", type=str, help="output text file for the covariance")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.cov_utils import export_cov

export_cov(args.rascalc_file, args.output_cov_file)