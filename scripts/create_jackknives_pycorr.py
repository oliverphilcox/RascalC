""" This script will assign a jackknife region to each input particle, in a way compatible with pycorr used in DESI data processing. Data is saved as a 5 column text file."""

import argparse

parser = argparse.ArgumentParser(description="This script will assign a jackknife region to each input particle, in a way compatible with pycorr used in DESI data processing. Data is saved as a 5 column text file.")
parser.add_argument("reference_file", type=str, help="reference (galaxies/data) text file name (rdz, i.e. RA, DEC, redshift columns; can have additional columns but they will not be used)")
parser.add_argument("input_file", type=str, help="input (randoms) text file name (xyzw, i.e. 3D Cartesian coordinates and weight columns)")
parser.add_argument("output_file", type=str, help="output (randoms) text file name (xyzwj, i.e. 3D Cartesian coordinates, weight and jackknife region number columns)")
parser.add_argument("n_jack", type=int, help="the number of jackknife regions to create")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pre_process.create_jackknives_pycorr import create_jackknives_pycorr_files

create_jackknives_pycorr_files(args.reference_file, args.input_file, args.output_file, args.n_jack)