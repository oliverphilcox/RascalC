### Script to smoothen a measured 2PCF in xi(r,mu) format.
### This computes all even multipoles xi_ell(r) up to a specified maximum ell, approximating the integral by a sum.
### Then r^2*xi_ell(r) are smoothened with Savitzky-Golay filter
### Finally, they are transformed back to xi(r,mu) format.

import argparse

parser = argparse.ArgumentParser(description="Script to smoothen a measured 2PCF in xi(r,mu) format. This computes all even multipoles xi_ell(r) up to a specified maximum ell, approximating the integral by a sum. Then r^2*xi_ell(r) are smoothened with Savitzky-Golay filter. Finally, they are transformed back to xi(r,mu) format.")
parser.add_argument("input_file", type=str, help="input xi text file")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("radial_window_length", type=int, help="width for the Savitzky-Golay filter")
parser.add_argument("radial_poly_order", type=int, help="polynomial order for the Savitzky-Golay filter")
parser.add_argument("output_file", type=str, help="output xi text file")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.xi.smoothen import smoothen_xi_files

smoothen_xi_files(args.input_file, args.max_l, args.radial_window_length, args.radial_pol_yorder, args.output_file)