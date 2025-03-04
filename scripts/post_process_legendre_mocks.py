## Script to post-process the single-field Legendre binned integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import sys

# PARAMETERS
if len(sys.argv) not in (6, 7, 8):
    print("Usage: python post_process_legendre_mocks.py {MOCK_COV_FILE} {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {OUTPUT_DIR} [{SKIP_R_BINS} [{SKIP_L}]]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process import post_process_legendre_mocks
        
mock_cov_file = str(sys.argv[1])
file_root = str(sys.argv[2])
n = int(sys.argv[3])
max_l = int(sys.argv[4])
outdir = str(sys.argv[5])
skip_r_bins = get_arg_safe(6, int, 0)
skip_l = get_arg_safe(7, int, 0)

post_process_legendre_mocks(mock_cov_file, file_root, n, max_l, outdir, skip_r_bins, skip_l)