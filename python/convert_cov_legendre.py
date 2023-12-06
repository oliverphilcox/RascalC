"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC Legendre mode results and saves full cov to text file, changing the indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys
from typing import Callable
from .convert_cov import load_cov, get_cov_header


def convert_cov_legendre(cov: np.ndarray[float], max_l: int) -> np.ndarray[float]:
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    n_bins = len(cov)
    if n_bins % n_l != 0: raise ValueError("Number of bins in the covariance must be divisible by the number of even multipoles")
    n_r_bins = n_bins // n_l
    cov = cov.reshape(n_r_bins, n_l, n_r_bins, n_l) # convert to 4D from 2D with [r, l] ordering for both rows and columns
    cov = cov.transpose(1, 0, 3, 2) # change orderng to [l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
    return cov

def load_cov_legendre(rascalc_results_file: str, max_l: int, print_function: Callable = print) -> np.ndarray[float]:
    return convert_cov_legendre(load_cov(rascalc_results_file, print_function), max_l)

def export_cov_legendre(rascalc_results_file: str, max_l: int, output_cov_file: str, print_function: Callable = print) -> None:
    np.savetxt(output_cov_file, load_cov_legendre(rascalc_results_file, max_l, print_function = print_function), header = get_cov_header(rascalc_results))

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python convert_cov_legendre.py {RASCALC_RESULTS_FILE} {MAX_L} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    max_l = int(sys.argv[2])
    output_cov_file = str(sys.argv[3])

    export_cov_legendre(rascalc_results, max_l, output_cov_file)
