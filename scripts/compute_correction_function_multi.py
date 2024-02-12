### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a two surveys. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import sys

# PARAMETERS
if len(sys.argv) not in (6, 9):
    print("Usage: python compute_correction_function_multi.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS_11} {RR_COUNTS_12} {RR_COUNTS_22}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.correction_function import compute_correction_function_multi

random_file = str(sys.argv[1])
random_file2 = str(sys.argv[2])
binfile = str(sys.argv[3])
outdir = str(sys.argv[4])
periodic = int(sys.argv[5])
RR_file = get_arg_safe(6)
RR_file12 = get_arg_safe(7)
RR_file2 = get_arg_safe(8)

compute_correction_function_multi(random_file, random_file2, binfile, outdir, periodic, RR_file, RR_file12, RR_file2)