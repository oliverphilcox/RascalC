### Script to fit a model (piecewise polynomial) to the 3PCF survey correction function, defined as the ratio between idealistic and true RRR pair counts for a single survey.

## NB: Input RRR counts should be normalized by the cube of the sum of random weights here. ## is this useful/necessary?
## NB: Assume mu is in [-1,1] limit here


import argparse

parser = argparse.ArgumentParser(description="Script to fit a model (piecewise polynomial) to the 3PCF survey correction function, defined as the ratio between idealistic and true RRR pair counts for a single survey.")
parser.add_argument("random_particle_file", type=str, help="file containing the random points' coordinates and weights")
parser.add_argument("r_bin_file", type=str, help="file containing the radial/separation bin boundaries in its rows")
parser.add_argument("output_dir", type=str, help="directory to write the resulting correction function coefficients")
parser.add_argument("periodic", type=bool, help="periodic boundary conditions flag")
parser.add_argument("RRR_counts_file", type=str, default=None, help="file containing the RRR counts in radial/separation and angular bins. necessary with aperiodic (realistic survey) geometry. not used with periodic geometry.")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.correction_function_3pcf import compute_3pcf_correction_function

compute_3pcf_correction_function(args.random_particle_file, args.r_bin_file, args.output_dir, args.periodic, args.RRR_counts_file)