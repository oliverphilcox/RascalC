"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC Legendre mode results and saves full cov to text file, changing the indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix"

import numpy as np
import sys
from convert_cov import load_cov


def convert_cov_legendre(cov: np.ndarray[float], n_r_bins: int):
    n_bins = len(cov)
    if n_bins % n_r_bins != 0: raise ValueError("Number of bins in the covariance must be divisible by the number of radial bins")
    n_l = n_bins // n_r_bins
    cov = cov.reshape(n_r_bins, n_l, n_r_bins, n_l) # convert to 4D from 2D with [r, l] ordering for both rows and columns
    cov = cov.transpose(1, 0, 3, 2) # change orderng to [l, r] for both rows and columns
    cov = cov.reshape(n_bins, n_bins) # convert back from 4D to 2D
    return cov

def load_cov_legendre(rascalc_results_file, n_r_bins, print_function = print):
    return convert_cov_legendre(load_cov(rascalc_results_file, print_function), n_r_bins)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python convert_cov_legendre.py {RASCALC_RESULTS_FILE} {N_R_BINS} {OUTPUT_COV_FILE}.")
        sys.exit(1)
    rascalc_results = str(sys.argv[1])
    n_r_bins = int(sys.argv[2])
    output_cov_file = str(sys.argv[3])

    np.savetxt(output_cov_file, load_cov_legendre(rascalc_results, n_r_bins))