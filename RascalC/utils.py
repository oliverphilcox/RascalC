# Contains some utility functions widely used in other scripts
# Not intended for execution from command line
import numpy as np
import os
from typing import Iterable


suffixes_tracer_all = ("", "2") # all supported tracer suffixes
indices_corr_all = ("11", "12", "22") # all supported 2PCF indices
suffixes_corr_all = ("", "12", "2") # all supported 2PCF suffixes
tracer1_corr = (0, 0, 1)
tracer2_corr = (0, 1, 1)


def blank_function(*args, **kwargs) -> None:
    """
    function that accepts anything and does nothing
    mostly intended for skipping optional printing
    """
    pass


def my_a2s(a, fmt='%.18e'):
    "custom array to string function"
    return ' '.join([fmt % e for e in a])


def transposed(A: np.ndarray):
    "swap last two (matrix) axes"
    return A.swapaxes(-2, -1)


def symmetrized(A: np.ndarray):
    "symmetrize a 2+D matrix over the last two axes"
    return 0.5 * (A + transposed(A))


def rmdir_if_exists_and_empty(dirname: str) -> None:
    "remove directory if it exists and is empty, otherwise do nothing"
    if os.path.isdir(dirname) and len(os.listdir(dirname)) <= 0:
        os.rmdir(dirname)


def format_skip_r_bins(skip_r_bins: int | tuple[int, int]) -> tuple[int, int]:
    """
    Format skip_r_bins into a tuple of two integers: number of bins to skip from the start and from the end.
    Single number results in skipping only from the start.
    """
    message = "`skip_r_bins` must be either an integer or a tuple of two integers"
    if isinstance(skip_r_bins, Iterable):
        if any(type(_) != int for _ in skip_r_bins): raise TypeError(message)
        if len(skip_r_bins) == 2: return tuple(skip_r_bins)
        if len(skip_r_bins) == 1: return skip_r_bins[0], 0
        if len(skip_r_bins) == 0: return 0, 0
        raise ValueError(message)
    if type(skip_r_bins) == int: return skip_r_bins, 0
    raise TypeError(message)