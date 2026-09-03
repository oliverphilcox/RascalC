"""Convenience script to normalize weights in an input (x,y,z,w) or (x,y,z,w,j) txt file and save the result in another .txt file for use with further scripts. (Michael Rashkovetskyi 2023).

    Parameters:
        INFILE = input ASCII file
        OUTFILE = output ASCII file

"""

import argparse

parser = argparse.ArgumentParser(description="Convenience script to normalize weights in an input (x,y,z,w) or (x,y,z,w,j) text file and save the result in another text file for use with further scripts.")
parser.add_argument("input_file", type=str)
parser.add_argument("output_file", type=str)
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pre_process.normalize_weights import normalize_weights_files

normalize_weights_files(args.input_file, args.output_file)