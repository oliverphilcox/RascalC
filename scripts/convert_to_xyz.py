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

import argparse

parser = argparse.ArgumentParser(description="Convenience script to convert an input (Ra,Dec,w) FITS or txt file to comoving (x,y,z) coordinates saved as a text file for use with the main C++ code. Output file format has (x,y,z,w) coordinates in Mpc/h units")
parser.add_argument("input_file", type=str, help="input text or FITS file name (text must be rdzw, i.e. RA, DEC, redshift and weight columns)")
parser.add_argument("output_file", type=str, help="output text file name (xyzw, i.e. 3D Cartesian coordinates and weight columns)")
parser.add_argument("Omega_m", type=float, default=0.31, nargs="?", help="relative matter density (default 0.31)")
parser.add_argument("Omega_K", type=float, default=0, nargs="?", help="relative curvature density (default 0)")
parser.add_argument("w_DE", type=float, default=-1, nargs="?", help="dark energy equation of state parameter (default -1)")
parser.add_argument("FKP_weights", type=str, default="0", nargs="?", help="whether to use FKP weights column (default False/0). also accepts format P0,NZ_name. only applies to (DESI) FITS files")
parser.add_argument("mask", type=int, default=0, nargs="?", help="sets bins that all must be set in STATUS for the particle to be selected (default 0). only applies to (DESI) FITS files")
parser.add_argument("use_weights", type=bool, default=True, nargs="?", help="whether to use WEIGHTS column. if not, set unit weights. default True/1")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pre_process.convert_to_xyz import convert_to_xyz_files
from RascalC.pre_process.utils import parse_FKP_arg

FKP_weights = parse_FKP_arg(args.FKP_weights) # FKP weights needs to be parsed additionally for backward compatibility (could be redesigned better)

convert_to_xyz_files(args.input_file, args.output_file, args.Omega_m, args.Omega_K, args.w_DE, FKP_weights, args.mask, args.use_weights)