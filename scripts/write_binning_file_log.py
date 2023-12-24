## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys

if len(sys.argv) < 5:
    print("Please specify input parameters in the form {N_LOG_BINS} {MIN_R} {MAX_R} {OUTPUT_FILE}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.write_binning_file import write_binning_file_log

nrbins = int(sys.argv[1])
r_min = float(sys.argv[2])
r_max = float(sys.argv[3])
out_file = str(sys.argv[4])

write_binning_file_log(out_file, r_min, r_max, nrbins)