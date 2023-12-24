"This reads a cosmodesi/pycorr .npy file and generates jackknife xi, weight, paircounts, and total paircounts text files for RascalC to use"

import numpy as np
import sys

## PARAMETERS
if len(sys.argv) not in (6, 7, 8, 9, 10, 11):
    print("Usage: python convert_xi_jack_from_pycorr.py {INPUT_NPY_FILE} {OUTPUT_XI_JACK_FILE} {OUTPUT_JACKKNIFE_WEIGHTS_FILE} {OUTPUT_JACKKNIFE_PAIRCOUNTS_FILE} {OUTPUT_BINNED_PAIRCOUNTS_FILE} [{R_STEP} [{N_MU} [{COUNTS_FACTOR} [{SPLIT_ABOVE} [{R_MAX}]]]]]. COUNTS_FACTOR=0 is a special value to use normalized counts.")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.pycorr_utils.jack import convert_jack_xi_weights_counts_from_pycorr_files

infile_name = str(sys.argv[1])
xi_jack_name = str(sys.argv[2])
jackweights_name = str(sys.argv[3])
jackpairs_name = str(sys.argv[4])
binpairs_name = str(sys.argv[5])
r_step = get_arg_safe(6, float, 1)
n_mu = get_arg_safe(7, int, None)
counts_factor = get_arg_safe(8, float, None) # basically number of randoms used for these counts, used to convert from total to 1 catalog count estimate
split_above = get_arg_safe(9, float, np.inf) # divide weighted RR counts by counts_factor**2 below this and by counts_factor above
r_max = get_arg_safe(10, float, np.inf) # if given, limit used r_bins at that

convert_jack_xi_weights_counts_from_pycorr_files(infile_name, xi_jack_name, jackweights_name, jackpairs_name, binpairs_name, n_mu, r_step, r_max, counts_factor, split_above)