## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys

if len(sys.argv) < 7:
    print("Please specify input parameters in the form {N_LOG_BINS} {N_LIN_BINS} {MIN_R} {CUTOFF_R} {MAX_R} {OUTPUT_FILE}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.write_binning_file import write_binning_file_hybrid

n_log_bins = int(sys.argv[1])
n_lin_bins = int(sys.argv[2])
r_min = float(sys.argv[3])
r_cut = float(sys.argv[4])
r_max = float(sys.argv[5])
out_file = str(sys.argv[6])

write_binning_file_hybrid(out_file, r_min, r_cut, r_max, n_log_bins, n_lin_bins)