## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys
import numpy as np
from utils import write_binning_file


def write_binning_file_log(out_file: str, r_min: float, r_cut: float, r_max: float, n_log_bins: int, n_lin_bins: int, print_function = print):
    print_function("Using hybrid binning: log binning up to r = %.1f, linear binning up to r = %.1f" %(r_cut,r_max))
    if r_min <= 0: raise ValueError('Minimum r must be positive to take logarithm')
    # Define radial bins
    r_edges = np.concatenate((np.geomspace(r_min, r_cut, n_log_bins+1)[:-1], np.linspace(r_cut, r_max, n_lin_bins+1)))
    write_binning_file(out_file, r_edges, print_function)

if __name__ == "__main__": # if invoked as a script
    if len(sys.argv)<7:
        print("Please specify input parameters in the form {N_LOG_BINS} {N_LIN_BINS} {MIN_R} {CUTOFF_R} {MAX_R} {OUTPUT_FILE}.")
        sys.exit(1)
    n_log_bins = int(sys.argv[1])
    n_lin_bins = int(sys.argv[2])
    r_min = float(sys.argv[3])
    r_cut = float(sys.argv[4])
    r_max = float(sys.argv[5])
    out_file = str(sys.argv[6])

    write_binning_file_log(out_file, r_min, r_cut, r_max, n_log_bins, n_lin_bins)