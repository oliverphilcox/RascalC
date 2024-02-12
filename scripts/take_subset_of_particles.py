"""This script takes an input txt file listing the (x,y,z,w,j) coordinates of random particles and selects N random particles from this, which it writes to file."""

import sys

if len(sys.argv) != 4:
    print("Usage: python take_subset_of_particles.py {INPUT_FILE} {OUTPUT_FILE} {N_PARTICLES}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.pre_process.take_subset_of_particles import take_subset_of_particles

infile_name = str(sys.argv[1])
outfile_name = str(sys.argv[2])
N = int(sys.argv[3])

take_subset_of_particles(infile_name, outfile_name, N)