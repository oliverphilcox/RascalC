## Simple convenience script to print shot noise rescaling from one or many RascalC files

import sys
import numpy as np
from .utils import blank_function


def get_shot_noise_rescaling(rascalc_filename: str) -> float | np.ndarray[float]:
    with np.load(rascalc_filename) as f:
        return f['shot_noise_rescaling']

def get_shot_noise_rescalings(rascalc_filenames: list[str], print_function = blank_function) -> dict:
    result = dict()
    for rascalc_filename in rascalc_filenames:
        result[rascalc_filename] = get_shot_noise_rescaling(rascalc_filename)
        print_function(f"{rascalc_filename}: shot-noise rescaling = {result[rascalc_filename]}")
    return result

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) < 2:
        print("Usage: python print_shot_noise_rescaling.py {RASCALC_RESULTS_1} [{RASCALC_RESULTS_2} ...]")
        sys.exit(1)
    
    get_shot_noise_rescalings(sys.argv[1:], print)