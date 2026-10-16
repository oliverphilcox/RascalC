### Script to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a two-tracer survey. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import argparse

parser = argparse.ArgumentParser(description="Script to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a two-tracer survey.")
parser.add_argument("random_particle_files", type=str, nargs=2, help="2 files containing the random points' coordinates and weights (one file for each tracer)")
parser.add_argument("r_bin_file", type=str, help="file containing the radial/separation bin boundaries in its rows")
parser.add_argument("output_dir", type=str, help="directory to write the resulting correction function coefficients")
parser.add_argument("periodic", type=bool, help="periodic boundary conditions flag")
parser.add_argument("RR_counts_files", type=str, default=[None]*3, nargs='*', help="3 files containing the RR counts in radial/separation and angular bins: first auto-counts for the first tracer, then the cross-counts between the tracers, and finally auto-counts for the second tracer. necessary with aperiodic (realistic survey) geometry. not used with periodic geometry.")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.correction_function import compute_correction_function_multi_from_files

compute_correction_function_multi_from_files(args.random_particle_files[0], args.random_particle_files[1], args.r_bin_file, args.output_dir, args.periodic, args.RR_counts_files[0], args.RR_counts_files[1], args.RR_counts_files[2])