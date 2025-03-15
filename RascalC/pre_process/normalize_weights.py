"""Convenience script to normalize weights in an input (x,y,z,w) or (x,y,z,w,j) txt file and save the result in another .txt file for use with further scripts. (Michael Rashkovetskyi 2023).

    Parameters:
        INFILE = input ASCII file
        OUTFILE = output ASCII file

"""

import numpy as np
from typing import Callable


def normalize_weights_files(input_file: str, output_file: str, print_function: Callable[[str], None] = print) -> None:
    print_function("Reading from file %s" % input_file)
    contents = np.loadtxt(input_file)
    print_function("Read %d particles" % contents.shape[0])
    if contents.ndim != 2: raise ValueError("Input file not read as 2D array")
    if contents.shape[1] < 4: raise ValueError("Not enough columns to get weights")

    print_function("Normalizing weights")
    contents[:, 3] /= np.sum(contents[:, 3]) # weights are the 4th column in r,d,z,w, x,y,z,w and x,y,z,w,j formats

    print_function("Saving to file %s" % output_file)
    np.savetxt(output_file, contents) # format?
