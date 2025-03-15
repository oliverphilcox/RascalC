## Simple functions to get shot noise rescaling from one or many RascalC files

import numpy as np
from .utils import blank_function
from typing import Callable


def get_shot_noise_rescaling(rascalc_filename: str) -> float | np.ndarray[float]:
    "Retrieve the shot noise rescaling value from RascalC Numpy (.npz) file."
    with np.load(rascalc_filename) as f:
        return f['shot_noise_rescaling']


def get_shot_noise_rescalings(rascalc_filenames: list[str], print_function: Callable[[str], None] = blank_function) -> dict:
    """
    Retrieve the shot noise rescaling value from multiple RascalC Numpy (.npz) file as a dictionary with filenames as keys.
    Optionally, use print_function to report the results.
    """
    result = dict()
    for rascalc_filename in rascalc_filenames:
        result[rascalc_filename] = get_shot_noise_rescaling(rascalc_filename)
        print_function(f"{rascalc_filename}: shot-noise rescaling = {result[rascalc_filename]}")
    return result