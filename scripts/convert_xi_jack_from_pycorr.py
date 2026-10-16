"This reads a cosmodesi/pycorr .npy file and generates jackknife xi, weight, paircounts, and total paircounts text files for RascalC to use"

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="This script reads cosmodesi/pycorr .npy file(s) and generates a binned pair counts text file for RascalC to use")
parser.add_argument("pycorr_file", type=str, help="pycorr .npy filename")
parser.add_argument("output_xi_jack_file", type=str, help="output text file for the jackknife 2PCF (xi)")
parser.add_argument("output_jack_weights_file", type=str, help="output text file for the jackknife region weights")
parser.add_argument("output_jack_RR_file", type=str, help="output text file for the jackknife binned RR pair counts")
parser.add_argument("output_RR_file", type=str, help="output text file for the total binned RR pair counts")
parser.add_argument("r_step", type=float, default=1, nargs='?', help="desired width of radial bins (default 1)")
parser.add_argument("n_mu_bins", type=int, default=0, nargs='?', help="number of angular (mu) bins (default 0, which means to preserve the original number)")
parser.add_argument("counts_factor", type=float, default=0, nargs='?', help="basically number of randoms used for these counts, used to convert from total to 1 catalog count estimate. 0 is a special value to use normalized counts, which is easiest and also the default")
parser.add_argument("split_above", type=float, default=np.inf, nargs='?', help="divide weighted RR counts by counts_factor squared below this and by counts_factor above. default is infinity")
parser.add_argument("r_max", type=float, default=np.inf, nargs='?', help="maximum radius (cutoff), the default is infinity (no cutoff)")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pycorr_utils.jack import convert_jack_xi_weights_counts_from_pycorr_files

# x or y should be a shortcut for x if x else y (or more explicitly x if x != 0 else y)
convert_jack_xi_weights_counts_from_pycorr_files(args.pycorr_file, args.output_xi_jack_file, args.output_jack_weights_file, args.output_jack_RR_file, args.output_RR_file, args.n_mu_bins or None, args.r_step, args.r_max, args.counts_factor or None, args.split_above)