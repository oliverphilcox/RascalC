"""Convenience script to perform redshift cut on an input (Ra,Dec,z,w) FITS or txt file and save the result in a .txt file for use with further scripts. (Michael Rashkovetskyi 2022).
Output file format has (Ra,Dec,z,w) coordinates

    Parameters:
        INFILE = input ASCII or FITS file
        OUTFILE = output .txt or .csv file specifier
        Z_MIN = Minimum redshift (inclusive, i.e. will require Z >= Z_MIN)
        Z_MAX = Maximum redshift (non-inclusive, i.e. will require Z < Z_MAX)
        ---OPTIONAL---
        USE_FKP_WEIGHTS = whether to use FKP weights column (default False/0; only applies to (DESI) FITS files)
        MASK = sets bins that all must be set in STATUS for the particle to be selected (default 0, only applies to (DESI) FITS files)
        USE_WEIGHTS = whether to use WEIGHTS column, if not, set unit weights (default True/1)

"""

import argparse

parser = argparse.ArgumentParser(description="Convenience script to perform redshift cut on an input (Ra,Dec,z,w) FITS or txt file and save the result in a text file for use with further scripts. Output file format has (Ra,Dec,z,w) coordinates")
parser.add_argument("input_file", type=str, help="input text or FITS file name (text must be rdzw, i.e. RA, DEC, redshift and weight columns)")
parser.add_argument("output_file", type=str, help="output text file name (rdzw, i.e. RA, DEC, redshift and weight columns)")
parser.add_argument("z_min", type=float)
parser.add_argument("z_max", type=float)
parser.add_argument("FKP_weights", type=str, default=0, nargs="?", help="whether to use FKP weights column (default False/0). also accepts format P0,NZ_name. only applies to (DESI) FITS files")
parser.add_argument("mask", type=int, default="0", nargs="?", help="sets bins that all must be set in STATUS for the particle to be selected (default 0). only applies to (DESI) FITS files")
parser.add_argument("use_weights", type=bool, default=True, nargs="?", help="whether to use WEIGHTS column. if not, set unit weights. default True/1")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pre_process.redshift_cut import redshift_cut_files
from RascalC.pre_process.utils import parse_FKP_arg

FKP_weights = parse_FKP_arg(args.FKP_weights) # FKP weights needs to be parsed additionally for backward compatibility (could be redesigned better)

redshift_cut_files(args.input_file, args.output_file, args.z_min, args.z_max, FKP_weights, args.mask, args.use_weights)