## Script to collect the raw covariance matrices from the output directory of the C++ code

import sys

# PARAMETERS
if len(sys.argv) != 2: # if too few
    print("Usage: python collect_raw_covariance_matrices.py {COVARIANCE_DIR}")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.raw_covariance_matrices import collect_raw_covariance_matrices

cov_dir = str(sys.argv[1])

collect_raw_covariance_matrices(cov_dir)