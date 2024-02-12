### Script to convert a measured 2PCF in xi(r,mu) format to Legendre multipoles, i.e. xi_ell(r).
### This computes all even multipoles up to a specified maximum ell, approximating the integral by a sum.
### The output form is a text file with the first column specifying the r-bin, the second giving xi_0(r), the third with xi_2(r) etc.

import sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python convert_xi_to_multipoles.py {INFILE} {MAX_L} {OUTFILE}")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.xi.convert_to_multipoles import convert_xi_to_multipoles_files

infile = str(sys.argv[1])
max_l = int(sys.argv[2])
outfile = str(sys.argv[3])

convert_xi_to_multipoles_files(infile, max_l, outfile)