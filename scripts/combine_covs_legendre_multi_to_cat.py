"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two full 2-tracer covs following NS/GCcomb procedure for 2 tracers in Legendre mode and convert them . Covariance of N(GC) and S(GC) 2PCF is neglected."

import sys

## PARAMETERS
if len(sys.argv) not in (13, 15, 17):
    print("Usage: python combine_covs_multi_to_cat.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1_11} {PYCORR_FILE2_11} {PYCORR_FILE1_12} {PYCORR_FILE2_12} {PYCORR_FILE1_22} {PYCORR_FILE2_22} {R_STEP} {MAX_L} {R_BINS_SKIP} {OUTPUT_COV_FILE} [{BIAS1} {BIAS2} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}]].")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.comb.combine_covs_legendre_multi_to_cat import combine_covs_legendre_multi_to_cat

rascalc_results1 = str(sys.argv[1])
rascalc_results2 = str(sys.argv[2])
pycorr_files1 = [str(sys.argv[i]) for i in (3, 5, 7)]
pycorr_files2 = [str(sys.argv[i]) for i in (4, 6, 8)]
r_step = float(sys.argv[9])
max_l = int(sys.argv[10])
skip_r_bins = int(sys.argv[11])
output_cov_file = str(sys.argv[12])
bias1 = get_arg_safe(13, float, 1)
bias2 = get_arg_safe(14, float, 1)
output_cov_file1 = get_arg_safe(15, str, None)
output_cov_file2 = get_arg_safe(16, str, None)

combine_covs_legendre_multi_to_cat(rascalc_results1, rascalc_results2, pycorr_files1, pycorr_files2, output_cov_file, max_l, r_step, skip_r_bins, bias1, bias2, output_cov_file1, output_cov_file2)