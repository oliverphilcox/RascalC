### Script to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a single survey. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import argparse

parser = argparse.ArgumentParser(description="Script to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a single survey.")
parser.add_argument("random_particle_file", type=str, help="file containing the random points' coordinates and weights")
parser.add_argument("r_bin_file", type=str, help="file containing the radial/separation bin boundaries in its rows")
parser.add_argument("output_dir", type=str, help="directory to write the resulting correction function coefficients")
parser.add_argument("periodic", type=bool, help="periodic boundary conditions flag")
parser.add_argument("RR_counts_file", type=str, default=None, nargs='?', help="file containing the RR counts in radial/separation and angular bins. necessary with aperiodic (realistic survey) geometry. not used with periodic geometry.")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.correction_function import compute_correction_function_from_file

compute_correction_function_from_file(args.random_particle_file, args.r_bin_file, args.output_dir, args.periodic, args.RR_counts_file)