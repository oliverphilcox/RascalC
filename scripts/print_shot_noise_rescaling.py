## Simple convenience script to print shot noise rescaling from one or many RascalC files

import sys

# PARAMETERS
if len(sys.argv) < 2:
    print("Usage: python print_shot_noise_rescaling.py {RASCALC_RESULTS_1} [{RASCALC_RESULTS_2} ...]")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.get_shot_noise_rescaling import get_shot_noise_rescalings

get_shot_noise_rescalings(sys.argv[1:], print_function = print)