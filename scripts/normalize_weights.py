"""Convenience script to normalize weights in an input (x,y,z,w) or (x,y,z,w,j) txt file and save the result in another .txt file for use with further scripts. (Michael Rashkovetskyi 2023).

    Parameters:
        INFILE = input ASCII file
        OUTFILE = output ASCII file

"""

import sys

# Check number of parameters
if len(sys.argv) != 3:
    print("Usage: python normalize_weights.py {INFILE} {OUTFILE}")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.pre_process.normalize_weights import normalize_weights_files

# Load file names
input_file = str(sys.argv[1])
output_file = str(sys.argv[2])

normalize_weights_files(input_file, output_file)