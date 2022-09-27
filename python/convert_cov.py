"Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of bias matrix"

import numpy as np
import sys

## PARAMETERS
if len(sys.argv) != 3:
    print("Usage: python convert_cov.py {RASCALC_RESULTS_FILE} {OUTPUT_COV_FILE}.")
    sys.exit()
rascalc_results = str(sys.argv[1])
output_cov_file = str(sys.argv[2])

with np.load(rascalc_results) as f:
    np.savetxt(output_cov_file, f['full_theory_covariance'])
    print(f"Max abs eigenvalue of bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
    # if the printed value is small the cov matrix should be safe to invert as is
