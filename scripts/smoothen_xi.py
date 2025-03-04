### Script to smoothen a measured 2PCF in xi(r,mu) format.
### This computes all even multipoles xi_ell(r) up to a specified maximum ell, approximating the integral by a sum.
### Then r^2*xi_ell(r) are smoothened with Savitzky-Golay filter
### Finally, they are transformed back to xi(r,mu) format.

import sys

## PARAMETERS
if len(sys.argv) != 6:
    print("Usage: python smoothen_xi.py {INFILE} {MAX_L} {RADIAL_WINDOW_LENGTH} {RADIAL_POLYORDER} {OUTFILE}")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.xi.smoothen import smoothen_xi_files

infile = str(sys.argv[1])
max_l = int(sys.argv[2])
window_length = int(sys.argv[3])
polyorder = int(sys.argv[4])
outfile = str(sys.argv[5])

smoothen_xi_files(infile, max_l, window_length, polyorder, outfile)