## Utility script to create the Corrfunc radial binning file for a given number of logarithmic and then linear radial bins in a certain range of comoving radii.

import argparse

parser = argparse.ArgumentParser(description="Utility script to create the Corrfunc radial binning file with logarithmic bins at smaller radii and linear bins at larger.")
parser.add_argument("n_log_r_bins", type=int, help="number of logarithmic radial bins")
parser.add_argument("n_lin_r_bins", type=int, help="number of linear radial bins")
parser.add_argument("r_min", type=float)
parser.add_argument("r_cut", type=float, help="boundary between the logarithmic and linear bins")
parser.add_argument("r_max", type=float)
parser.add_argument("output_file", type=str)
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.write_binning_file import write_binning_file_hybrid

write_binning_file_hybrid(args.output_file, args.r_min, args.r_cut, args.r_max, args.n_log_r_bins, args.n_lin_r_bins)