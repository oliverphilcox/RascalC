## Script to post-process the single-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import sys

# PARAMETERS
if len(sys.argv) not in (5, 6, 7):
    print("Usage: python post_process_default.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process.default import post_process_default

file_root = str(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
outdir = str(sys.argv[4])
alpha = get_arg_safe(5, float, 1)
skip_r_bins = get_arg_safe(6, int, 0)

post_process_default(file_root, n, m, outdir, alpha, skip_r_bins)