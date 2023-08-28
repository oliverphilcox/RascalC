"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC 2-tracer Legendre mode results and saves full cov to text file, changing the indexing from [t, r, l] to [t, l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python convert_cov_legendre_multi.py {RASCALC_RESULTS_FILE} {N_R_BINS} {OUTPUT_COV_FILE}.")
    sys.exit(1)
rascalc_results = str(sys.argv[1])
n = int(sys.argv[2])
output_cov_file = str(sys.argv[3])

with np.load(rascalc_results) as f:
    cov = f['full_theory_covariance']
    print(f"Max abs eigenvalue of bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
    # if the printed value is small the cov matrix should be safe to invert as is

n_bins = len(cov)
assert n_bins % (3*n) == 0, "Number of bins mismatch"
n_l = n_bins // (3*n)
cov = cov.reshape(3, n, n_l, 3, n, n_l) # convert to 6D from 2D with [t, r, l] ordering for both rows and columns
cov = cov.transpose(0, 2, 1, 3, 5, 4) # change orderng to [t, l, r] for both rows and columns
cov = cov.reshape(n_bins, n_bins) # convert back from 6D to 2D
np.savetxt(output_cov_file, cov)