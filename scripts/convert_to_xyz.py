"""Convenience script to convert an input (Ra,Dec,w) FITS or txt file to comoving (x,y,z) coordinates saved as a .txt file for use with the main C++ code. (Oliver Philcox 2018 with modifications by Michael Rashkovetskyi 2022, using Daniel Eisenstein's 2015 WCDM Coordinate Converter).
Output file format has (x,y,z,w) coordinates in Mpc/h units 

    Parameters:
        INFILE = input ASCII or FITS file
        OUTFILE = output .txt or .csv file specifier
        ---OPTIONAL---
        OMEGA_M = Matter density (default 0.31)
        OMEGA_K = Curvature density (default 0.0)
        W_DARK_ENERGY = Dark Energy equation of state parameter (default -1.)
        ---FURTHER OPTIONAL---
        USE_FKP_WEIGHTS = whether to use FKP weights column (default False/0; only applies to (DESI) FITS files)
        MASK = sets bins that all must be set in STATUS for the particle to be selected (default 0, only applies to (DESI) FITS files)
        USE_WEIGHTS = whether to use WEIGHTS column, if not, set unit weights (default True/1)

"""

import sys

if len(sys.argv) not in (3, 4, 5, 6, 7, 8, 9):
    print("Usage: python convert_to_xyz.py {INFILE} {OUTFILE} [{OMEGA_M} {OMEGA_K} {W_DARK_ENERGY} [{USE_FKP_WEIGHTS or P0,NZ_name} [{MASK} [{USE_WEIGHTS}]]]]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.pre_process.convert_to_xyz import convert_to_xyz_files
from RascalC.pre_process.utils import parse_FKP_arg, my_str_to_bool
        
# Load file names
input_file = str(sys.argv[1])
output_file = str(sys.argv[2])

# Read in optional cosmology parameters
Omega_m = get_arg_safe(3, float, 0.31)
Omega_k = get_arg_safe(4, float, 0)
w_dark_energy = get_arg_safe(5, float, -1)
# defaults from the BOSS DR12 2016 clustering paper assuming LCDM

# The next only applies to (DESI) FITS files
FKP_weights = parse_FKP_arg(get_arg_safe(6, str, "False")) # determine whether to use FKP weights
mask = get_arg_safe(7, int, 0) # default is 0 - no mask filtering
use_weights = my_str_to_bool(get_arg_safe(8, str, "True")) # use weights by default

convert_to_xyz_files(input_file, output_file, Omega_m, Omega_k, w_dark_energy, FKP_weights, mask, use_weights)