"This reads a .npy file of RascalC Legendre results (or txt covariance converted previously) and a triplet of cosmodesi/pycorr .npy files to produce a covariance for a catalog of these two tracers concatenated."

import sys

## PARAMETERS
if len(sys.argv) not in (9, 11):
    print("Usage: python combine_cov_multi_to_cat.py {RASCALC_RESULTS} {PYCORR_FILE_11} {PYCORR_FILE_12} {PYCORR_FILE_22} {R_STEP} {MAX_L} {SKIP_R_BINS} {OUTPUT_COV_FILE} [{BIAS1} {BIAS2}].")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.comb.convert_cov_legendre_multi_to_cat import convert_cov_legendre_multi_to_cat

rascalc_results = str(sys.argv[1])
pycorr_files = [str(sys.argv[i]) for i in range(2, 5)]
r_step = float(sys.argv[5])
max_l = int(sys.argv[6])
skip_r_bins = int(sys.argv[7])
output_cov_file = str(sys.argv[8])
from .utils import get_arg_safe
bias1 = get_arg_safe(9, float, 1)
bias2 = get_arg_safe(10, float, 1)

convert_cov_legendre_multi_to_cat(rascalc_results, pycorr_files, output_cov_file, max_l, r_step, skip_r_bins, bias1, bias2)