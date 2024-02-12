### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a single survey. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import sys

# PARAMETERS
if len(sys.argv) not in (5, 6):
    print("Usage: python compute_correction_function.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS}] ")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.correction_function import compute_correction_function

random_file = str(sys.argv[1])
binfile = str(sys.argv[2])
outdir = str(sys.argv[3])
periodic = int(sys.argv[4])
RR_file = get_arg_safe(5)

compute_correction_function(random_file, binfile, outdir, periodic, RR_file)