""" This function will assign a jackknife region to each input particle, in a way compatible with pycorr used in DESI data processing. Data is saved as a 5 column text file."""

import sys

## PARAMETERS
if len(sys.argv) != 5:
    print("Usage: python create_jackknives_pycorr.py {REF_RDZ_FILE} {INPUT_XYZW_FILE} {OUTPUT_XYZWJ_FILE} {NJACK}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.pre_process.create_jackknives_pycorr import create_jackknives_pycorr_files

reffile_name = str(sys.argv[1])
infile_name = str(sys.argv[2])
outfile_name = str(sys.argv[3])
njack = int(sys.argv[4])

create_jackknives_pycorr_files(reffile_name, infile_name, outfile_name, njack)