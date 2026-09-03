### Script to convert a measured 2PCF in xi(r,mu) format to Legendre multipoles, i.e. xi_ell(r).
### This computes all even multipoles up to a specified maximum ell, approximating the integral by a sum.
### The output form is a text file with the first column specifying the r-bin, the second giving xi_0(r), the third with xi_2(r) etc.

import argparse

parser = argparse.ArgumentParser(description="Script to convert a measured 2PCF in xi(r,mu) format to Legendre multipoles, i.e. xi_ell(r). This computes all even multipoles up to a specified maximum ell, approximating the integral by a sum. The output form is a text file with the first column specifying the r-bin, the second giving xi_0(r), the third with xi_2(r) etc.")
parser.add_argument("input_file", type=str, help="input xi text file")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_file", type=str, help="output xi text file")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.xi.convert_to_multipoles import convert_xi_to_multipoles_files

convert_xi_to_multipoles_files(args.input_file, args.max_l, args.output_file)