## Utility script to create the Corrfunc radial binning file for a given number of logarithmic radial bins in a certain range of comoving radii.

import argparse

parser = argparse.ArgumentParser(description="Utility script to create the Corrfunc radial binning file with logarithmic bins.")
parser.add_argument("n_r_bins", type=int, help="number of radial bins")
parser.add_argument("r_min", type=float)
parser.add_argument("r_max", type=float)
parser.add_argument("output_file", type=str)
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.write_binning_file import write_binning_file_log

write_binning_file_log(args.output_file, args.r_min, args.r_max, args.n_r_bins)