"This reads a cosmodesi/pycorr .npy file and generates binned pair counts text file for RascalC to use"

import sys
import numpy as np

## PARAMETERS
if len(sys.argv) not in (3, 4, 5, 6, 7, 8):
    print("Usage: python convert_counts_from_pycorr.py {INPUT_NPY_FILE} {OUTPUT_PAIRCOUNTS_TEXT_FILE} [{R_STEP} [{N_MU} [{COUNTS_FACTOR} [{SPLIT_ABOVE}] [{R_MAX}]]]]]. COUNTS_FACTOR=0 is a special value to use normalized counts.")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.pycorr_utils.counts import convert_counts_from_pycorr_files

infile_name = str(sys.argv[1])
outfile_name = str(sys.argv[2])
r_step = get_arg_safe(3, float, 1)
n_mu = get_arg_safe(4, int, None)
counts_factor = get_arg_safe(5, float, None) # basically number of randoms used for these counts, used to convert from total to 1 catalog count estimate
split_above = get_arg_safe(6, float, np.inf) # divide weighted RR counts by counts_factor**2 below this and by counts_factor above
r_max = get_arg_safe(7, float, np.inf) # if given, limit used r_bins at that

if n_mu == 0: n_mu = None
if counts_factor == 0: counts_factor = None
convert_counts_from_pycorr_files(infile_name, outfile_name, n_mu, r_step, r_max, counts_factor, split_above)