"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC Legendre mode results and saves full cov to text file, changing the indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python convert_cov_legendre.py {RASCALC_RESULTS_FILE} {N_R_BINS} {OUTPUT_COV_FILE}.")
    sys.exit()
rascalc_results = str(sys.argv[1])
n = int(sys.argv[2])
output_cov_file = str(sys.argv[3])

with np.load(rascalc_results) as f:
    cov = f['full_theory_covariance']
    print(f"Max abs eigenvalue of bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
    # if the printed value is small the cov matrix should be safe to invert as is

n_bins = len(cov)
assert n_bins % n == 0, "Number of bins mismatch"
n_l = n_bins // n
cov = cov.reshape(n, n_l, n, n_l) # convert to 4D from 2D with [r, l] ordering for both rows and columns
cov = cov.transpose(1, 0, 3, 2) # change orderng to [l, r] for both rows and columns
cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
np.savetxt(output_cov_file, cov)