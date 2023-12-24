## Script to catenate (take a subset of subsamples or resample) raw covariance matrix files produced by collect_raw_covariance_matrices.py
## More specifically, copy part of partial results to other directory and recompute totals by averaging
## Determines single-field vs multi-field and jackknife automatically
## Do not use if subsamples have different numbers of pairs/triples/quadruplets

import sys

# PARAMETERS
if len(sys.argv)<6: # if too few
    print("Usage: python cat_raw_covariance_matrices.py {N_R_BINS} {mN_MU_BINS/lMAX_L} {COVARIANCE_INPUT_DIR1} {N_SUBSAMPLES_TO_USE1} [{COVARIANCE_INPUT_DIR2} {N_SUBSAMPLES_TO_USE2} ...] [{COLLAPSE_FACTOR}] {COVARIANCE_OUTPUT_DIR}")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.raw_covariance_matrices import cat_raw_covariance_matrices

n = int(sys.argv[1])
mstr = str(sys.argv[2])
input_roots = [str(s) for s in sys.argv[3:-1:2]]
ns_samples = [int(s) for s in sys.argv[4:-1:2]]
output_root = str(sys.argv[-1])
collapse_factor = int(input_roots.pop()) if len(input_roots) > len(ns_samples) else 1 # recover the collapse factor if present

cat_raw_covariance_matrices(n, mstr, input_roots, ns_samples, output_root, collapse_factor)