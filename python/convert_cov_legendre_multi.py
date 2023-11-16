"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC 2-tracer Legendre mode results and saves full cov to text file, changing the indexing from [t, r, l] to [t, l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys
from convert_cov import load_cov


def convert_cov_legendre_multi(cov: np.ndarray[float], n_r_bins: int):
    n_bins = len(cov)
    if n_bins % (3 * n_r_bins) != 0: raise ValueError("Number of bins in the covariance must be divisible by thrice the number of radial bins")
    n_l = n_bins // (3 * n_r_bins)
    cov = cov.reshape(3, n_r_bins, n_l, 3, n_r_bins, n_l) # convert to 6D from 2D with [t, r, l] ordering for both rows and columns
    cov = cov.transpose(0, 2, 1, 3, 5, 4) # change orderng to [t, l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 6D to 2D
    return cov

def load_cov_legendre_multi(rascalc_results_file, n_r_bins, print_function = print):
    return convert_cov_legendre_multi(load_cov(rascalc_results_file, print_function), n_r_bins)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python convert_cov_legendre_multi.py {RASCALC_RESULTS_FILE} {N_R_BINS} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    n_r_bins = int(sys.argv[2])
    output_cov_file = str(sys.argv[3])

    np.savetxt(output_cov_file, load_cov_legendre_multi(rascalc_results, n_r_bins))