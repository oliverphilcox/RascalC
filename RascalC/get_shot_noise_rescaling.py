## Simple functions to get shot noise rescaling from one or many RascalC files

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