## Script to post-process the single-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import sys

# PARAMETERS
if len(sys.argv) not in (7, 8, 9):
    print("Usage: python post_process_legendre_mix_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {MAX_L} {OUTPUT_DIR} [{SKIP_R_BINS} [{SKIP_L}]]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.post_process.legendre_mix_jackknife import post_process_legendre_mix_jackknife
        
jackknife_file = str(sys.argv[1])
weight_dir = str(sys.argv[2])
file_root = str(sys.argv[3])
m = int(sys.argv[4])
max_l = int(sys.argv[5])
outdir = str(sys.argv[6])
skip_r_bins = get_arg_safe(7, int, 0)
skip_l = get_arg_safe(8, int, 0)

post_process_legendre_mix_jackknife(jackknife_file, weight_dir, file_root, m, max_l, outdir, skip_r_bins, skip_l)