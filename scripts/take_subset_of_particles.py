"""This script takes an input txt file listing the (x,y,z,w,j) coordinates of random particles and selects N random particles from this, which it writes to file."""

import argparse

parser = argparse.ArgumentParser(description="This script takes an input txt file listing the (x,y,z,w,j) coordinates of random particles, selects N random particles and writes them to the output file.")
parser.add_argument("input_file", type=str)
parser.add_argument("output_file", type=str)
parser.add_argument("N", type=int, help="number of particles to randomly select")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pre_process.take_subset_of_particles import take_subset_of_particles

take_subset_of_particles(args.input_file, args.output_file, args.N)