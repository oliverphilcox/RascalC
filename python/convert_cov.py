"Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of bias matrix"

import numpy as np
import sys


def load_cov(rascalc_results_file: str, print_function = print) -> np.ndarray[float]:
    with np.load(rascalc_results_file) as f:
        print_function(f"Max abs eigenvalue of bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
        # if the printed value is small the cov matrix should be safe to invert as is
        return f['full_theory_covariance']

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 3:
        print("Usage: python convert_cov.py {RASCALC_RESULTS_FILE} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    output_cov_file = str(sys.argv[2])

    np.savetxt(output_cov_file, load_cov(rascalc_results))