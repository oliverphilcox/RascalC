## Script to collect the raw covariance matrices from the output directory of the C++ code

import argparse

parser = argparse.ArgumentParser(description="Script to collect the raw covariance matrices from the output directory of the C++ code. It converts multiple text files to a single binary (Numpy) file.")
parser.add_argument("covariance_dir", type=str, help="directory containing the covariance matrix subdirectory")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.raw_covariance_matrices import collect_raw_covariance_matrices

collect_raw_covariance_matrices(args.covariance_dir)