"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC 2-tracer Legendre mode results and saves full cov to text file, changing the indexing from [t, r, l] to [t, l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys
from .convert_cov import load_cov, get_cov_header


def convert_cov_legendre_multi(cov: np.ndarray[float], max_l: int):
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    n_bins = len(cov)
    if n_bins % (3 * n_l) != 0: raise ValueError("Number of bins in the covariance must be divisible by thrice the number of even multipoles")
    n_r_bins = n_bins // (3 * n_l)
    cov = cov.reshape(3, n_r_bins, n_l, 3, n_r_bins, n_l) # convert to 6D from 2D with [t, r, l] ordering for both rows and columns
    cov = cov.transpose(0, 2, 1, 3, 5, 4) # change orderng to [t, l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 6D to 2D
    return cov

def load_cov_legendre_multi(rascalc_results_file, max_l, print_function = print):
    return convert_cov_legendre_multi(load_cov(rascalc_results_file, print_function), max_l)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python convert_cov_legendre_multi.py {RASCALC_RESULTS_FILE} {MAX_L} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    max_l = int(sys.argv[2])
    output_cov_file = str(sys.argv[3])

    np.savetxt(output_cov_file, load_cov_legendre_multi(rascalc_results, max_l), header = get_cov_header(rascalc_results))
