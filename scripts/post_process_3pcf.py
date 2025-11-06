## Script to post-process the single-field 3PCF Legendre binned integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import sys

# PARAMETERS
if len(sys.argv) not in (6, 7):
    print("Usage: python post_process_3pcf.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process_3pcf import post_process_3pcf
        
file_root = str(sys.argv[1])
n = int(sys.argv[2])
max_l = int(sys.argv[3])
n_samples = int(sys.argv[4])
outdir = str(sys.argv[5])
alpha = get_arg_safe(6, float, 1)

post_process_3pcf(file_root, n, max_l, n_samples, outdir, alpha)