## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys
import numpy as np
from .utils import write_binning_file


def write_binning_file_linear(out_file: str, r_min: float, r_max: float, nrbins: int, print_function = print):
    print_function("Using LINEAR binning")
    # Define radial bins
    r_edges = np.linspace(r_min, r_max, nrbins+1)
    if r_edges[0] <= 1e-4: r_edges[0] = 1e-4 # exclude very small separations for stability
    write_binning_file(out_file, r_edges, print_function)

if __name__ == "__main__": # if invoked as a script
    if len(sys.argv)<5:
        print("Please specify input parameters in the form {N_BINS} {MIN_R} {MAX_R} {OUTPUT_FILE}.")
        sys.exit(1)
    nrbins = int(sys.argv[1])
    r_min = float(sys.argv[2])
    r_max = float(sys.argv[3])
    out_file = str(sys.argv[4])

    write_binning_file_linear(out_file, r_min, r_max, nrbins)