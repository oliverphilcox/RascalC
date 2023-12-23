"Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of bias matrix"

import numpy as np
import sys
from typing import Callable
from .print_shot_noise_rescaling import get_shot_noise_rescaling


def get_cov_header(rascalc_results_file: str) -> str:
    return "shot_noise_rescaling = " + str(get_shot_noise_rescaling(rascalc_results_file))

def load_cov(rascalc_results_file: str, print_function: Callable = print) -> np.ndarray[float]:
    with np.load(rascalc_results_file) as f:
        print_function(f"Max abs eigenvalue of bias correction matrix is {np.max(np.abs(np.linalg.eigvals(f['full_theory_D_matrix']))):.2e}")
        # if the printed value is small the cov matrix should be safe to invert as is
        return f['full_theory_covariance']

def export_cov(rascalc_results_file: str, output_cov_file: str, print_function: Callable = print) -> None:
    np.savetxt(output_cov_file, load_cov(rascalc_results_file, print_function = print_function), header = get_cov_header(rascalc_results_file))

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 3:
        print("Usage: python convert_cov.py {RASCALC_RESULTS_FILE} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    output_cov_file = str(sys.argv[2])

    export_cov(rascalc_results, output_cov_file)