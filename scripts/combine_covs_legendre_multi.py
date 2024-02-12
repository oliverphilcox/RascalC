"This reads two sets of RascalC results and two triplets of cosmodesi/pycorr .npy files to combine two full 2-tracer covs following NS/GCcomb procedure for 2 tracers in Legendre mode. Covariance of N(GC) and S(GC) 2PCF is neglected."

import sys

## PARAMETERS
if len(sys.argv) not in (13, 15):
    print("Usage: python combine_covs_legendre_multi.py {RASCALC_RESULTS1} {RASCALC_RESULTS2} {PYCORR_FILE1_11} {PYCORR_FILE2_11} {PYCORR_FILE1_12} {PYCORR_FILE2_12} {PYCORR_FILE1_22} {PYCORR_FILE2_22} {N_R_BINS} {MAX_L} {SKIP_R_BINS} {OUTPUT_COV_FILE} [{OUTPUT_COV_FILE1} {OUTPUT_COV_FILE2}].")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.comb.combine_covs_legendre_multi import combine_covs_legendre_multi

rascalc_results1 = str(sys.argv[1])
rascalc_results2 = str(sys.argv[2])
pycorr_files1 = [str(sys.argv[i]) for i in (3, 5, 7)]
pycorr_files2 = [str(sys.argv[i]) for i in (4, 6, 8)]
r_step = float(sys.argv[9])
max_l = int(sys.argv[10])
skip_r_bins = int(sys.argv[11])
output_cov_file = str(sys.argv[12])
output_cov_file1 = get_arg_safe(13, str, None)
output_cov_file2 = get_arg_safe(14, str, None)

combine_covs_legendre_multi(rascalc_results1, rascalc_results2, pycorr_files1, pycorr_files2, output_cov_file, max_l, r_step, skip_r_bins, output_cov_file1, output_cov_file2)