"This reads two sets of RascalC results and two cosmodesi/pycorr .npy files to combine two covs following NS/GCcomb procedure in Legendre mode. Covariance of N and S 2PCF is neglected."

import sys

## PARAMETERS
if len(sys.argv) not in (9, 11):
    print("Usage: python combine_covs_legendre.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1} {PYCORR_FILE2} {R_STEP} {MAX_L} {SKIP_R_BINS} {OUTPUT_COV_FILE} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}].")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.comb.combine_covs_legendre import combine_covs_legendre

rascalc_results1 = str(sys.argv[1])
rascalc_results2 = str(sys.argv[2])
pycorr_file1 = str(sys.argv[3])
pycorr_file2 = str(sys.argv[4])
r_step = float(sys.argv[5])
max_l = int(sys.argv[6])
skip_r_bins = int(sys.argv[7])
output_cov_file = str(sys.argv[8])
output_cov_file1 = get_arg_safe(9, str, None)
output_cov_file2 = get_arg_safe(10, str, None)

combine_covs_legendre(rascalc_results1, rascalc_results2, pycorr_file1, pycorr_file2, output_cov_file, max_l, r_step, skip_r_bins, output_cov_file1, output_cov_file2)