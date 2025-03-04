## Script to post-process the multi-field integrals computed by the C++ code. This computes two shot-noise rescaling parameters, alphas, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import sys

# PARAMETERS
if len(sys.argv) not in (7, 8):
    print("Usage: python post_process_default_mocks_multi.py {MOCK_COV_FILE} {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SKIP_R_BINS}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process import post_process_default_mocks_multi
        
mock_cov_file = str(sys.argv[1])
file_root = str(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])
n_samples = int(sys.argv[5])
outdir = str(sys.argv[6])
skip_r_bins = get_arg_safe(7, int, 0) # convert from radial to total number of bins right away

post_process_default_mocks_multi(mock_cov_file, file_root, n, m, n_samples, outdir, skip_r_bins)
