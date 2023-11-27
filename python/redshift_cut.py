"""Convenience script to perform redshift cut on an input (Ra,Dec,z,w) FITS or txt file and save the result in a .txt file for use with further scripts. (Michael Rashkovetskyi 2022).
Output file format has (Ra,Dec,z,w) coordinates

    Parameters:
        INFILE = input ASCII or FITS file
        OUTFILE = output .txt or .csv file specifier
        Z_MIN = Minimum redshift (inclusive, i.e. will require Z >= Z_MIN)
        Z_MAX = Maximum redshift (non-inclusive, i.e. will require Z < Z_MAX)
        ---OPTIONAL---
        USE_FKP_WEIGHTS = whether to use FKP weights column (default False/0; only applies to (DESI) FITS files)
        MASK = sets bins that all must be set in STATUS for the particle to be selected (default 0, only applies to (DESI) FITS files)
        USE_WEIGHTS = whether to use WEIGHTS column, if not, set unit weights (default True/1)

"""

import sys
import numpy as np
from .utils import read_particles_fits_file


def redshift_cut_files(input_file: str, output_file: str, z_min: float, z_max: float, FKP_weights: bool | tuple[float, str] = False, mask: int = 0, use_weights: bool = True, print_function = print):
    # Load in data:
    print_function("Reading input file %s in Ra,Dec,z coordinates\n"%input_file)
    if input_file.endswith(".fits"):
        particles = read_particles_fits_file(input_file, FKP_weights, mask, use_weights)
    else:
        particles = np.loadtxt(input_file)
    # in either case, first index (rows) are particle numbers, second index (columns) are different properties

    print_function("Performing redshift cut")
    all_z = particles[:, 2] # assume the redshift is the third column
    filt = np.logical_and(z_min <= all_z, all_z < z_max)

    print_function("Writing to file %s:"%output_file)
    np.savetxt(output_file, particles[filt])
    print_function("Output positions (of length %d) written succesfully!" % len(particles))

if __name__ == "__main__": # if invoked as a script
    # Check number of parameters
    if len(sys.argv) not in (5, 6, 7, 8):
        print("Usage: python redshift_cut.py {INFILE} {OUTFILE} {Z_MIN} {Z_MAX} [{USE_FKP_WEIGHTS or P0,NZ_name} [{MASK} [{USE_WEIGHTS}]]]")
        sys.exit(1)
            
    # Load file names
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])
    print("\nUsing input file %s in Ra,Dec,z coordinates\n"%input_file)

    # Load min and max redshifts
    z_min = float(sys.argv[3])
    z_max = float(sys.argv[4])

    from .utils import get_arg_safe, parse_FKP_arg, my_str_to_bool
    # The next only applies to (DESI) FITS files
    FKP_weights = parse_FKP_arg(get_arg_safe(5, str, "False")) # determine whether to use FKP weights
    mask = get_arg_safe(6, int, 0) # default is 0 - no mask filtering
    use_weights = my_str_to_bool(get_arg_safe(7, str, "True")) # use weights by default

    redshift_cut_files(input_file, output_file, z_min, z_max, FKP_weights, mask, use_weights)