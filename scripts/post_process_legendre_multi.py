## Script to post-process the multi-field Legendre binned integrals computed by the C++ code, given a shot-noise rescaling parameter alpha.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import sys

# PARAMETERS
if len(sys.argv) not in (5, 7, 8, 9):
    print("Usage: python post_process_legendre_multi.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process import post_process_legendre_multi

file_root = str(sys.argv[1])
n = int(sys.argv[2])
max_l = int(sys.argv[3])
outdir = str(sys.argv[4])
alpha_1 = get_arg_safe(5, float, 1)
alpha_2 = get_arg_safe(6, float, 1)
skip_r_bins = get_arg_safe(7, int, 0)
skip_l = get_arg_safe(8, int, 0)

post_process_legendre_multi(file_root, n, max_l, outdir, alpha_1, alpha_2, skip_r_bins, skip_l)
