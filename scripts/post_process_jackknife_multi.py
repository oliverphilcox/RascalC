## Script to post-process the multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameters, alpha_i, from data derived covariance matrices.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices.
## The effective number of samples, N_eff, is also computed.

import sys

# PARAMETERS
if len(sys.argv) not in (8, 9):
    print("Usage: python post_process_jackknife_multi.py {XI_JACKKNIFE_FILE_11} {XI_JACKKNIFE_FILE_12} {XI_JACKKNIFE_FILE_22} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {OUTPUT_DIR} [{SKIP_R_BINS}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process.jackknife_multi import post_process_jackknife_multi

jackknife_file_11 = str(sys.argv[1])
jackknife_file_12 = str(sys.argv[2])
jackknife_file_22 = str(sys.argv[3])
weight_dir = str(sys.argv[4])
file_root = str(sys.argv[5])
m = int(sys.argv[6])
outdir = str(sys.argv[7])
skip_r_bins = get_arg_safe(8, int, 0)

post_process_jackknife_multi(jackknife_file_11, jackknife_file_12, jackknife_file_22, weight_dir, file_root, m, outdir, skip_r_bins)
