## Simple convenience script to print shot noise rescaling from one or many RascalC files

import argparse

parser = argparse.ArgumentParser(description="Simple convenience script to print shot noise rescaling from one or many RascalC files")
parser.add_argument("rascalc_results", type=str, nargs="+", help="RascalC .npy filename(s)")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.get_shot_noise_rescaling import get_shot_noise_rescalings

get_shot_noise_rescalings(args.rascalc_results, print_function = print)