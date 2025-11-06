### Script to fit a model (piecewise polynomial) to the 3PCF survey correction function, defined as the ratio between idealistic and true RRR pair counts for a single survey.

## NB: Input RRR counts should be normalized by the cube of the sum of random weights here. ## is this useful/necessary?
## NB: Assume mu is in [-1,1] limit here


import sys

# PARAMETERS
if len(sys.argv) not in (5, 6):
    print("Usage: python compute_3pcf_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RRR_COUNTS}]")
    sys.exit(1)

gal_file = str(sys.argv[1])
binfile = str(sys.argv[2])
outdir = str(sys.argv[3])
periodic = int(sys.argv[4])

from utils import get_arg_safe, adjust_path
RRR_file = get_arg_safe(5)
adjust_path()
from RascalC.correction_function_3pcf import compute_3pcf_correction_function

compute_3pcf_correction_function(gal_file, binfile, outdir, periodic, RRR_file)