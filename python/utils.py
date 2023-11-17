# Contains some utility functions widely used in other scripts
# Not intended for execution from command line
import sys, os
import numpy as np
from warnings import warn
from astropy.io import fits

def get_arg_safe(index: int, type = str, default: object = None) -> object:
    # get argument by index from sys.argv and convert it to the requested type if there are enough elements there
    # otherwise return the default value
    return type(sys.argv[index]) if len(sys.argv) > index else default

def blank_function(*args, **kwargs) -> None:
    # function that accepts anything and does nothing
    # mostly intended for skipping optional printing
    pass

def my_a2s(a, fmt='%.18e'):
    # custom array to string function
    return ' '.join([fmt % e for e in a])

def my_str_to_bool(s: str) -> bool:
    # naive conversion to bool for all non-empty strings is True, and one can't give an empty string as a command line argument, so need to make it more explicit
    return s not in ("0", "false")

def symmetrized(A):
    # symmetrize a 2D matrix
    return 0.5 * (A + A.T)

def parse_FKP_arg(FKP_weights: str) -> bool | tuple[float, str]:
    if not my_str_to_bool(FKP_weights): return False
    # determine if it actually has P0,NZ_name format. Such strings should convert to True
    arg_FKP_split = FKP_weights.split(",")
    if len(arg_FKP_split) == 2:
        return (float(arg_FKP_split[0]), arg_FKP_split[1])
    if len(arg_FKP_split) == 1: return True
    raise ValueError("FKP parameter matched neither USE_FKP_WEIGHTS (true/false in any register or 0/1) nor P0,NZ_name (float and string without space).")

def read_particles_fits_file(input_file: str, FKP_weights: bool | (float, str) = False, mask: int = 0, use_weights: bool = True):
    # Read FITS file with particles. Can apply mask filtering and compute FKP weights in different ways. Works for DESI setups
    filt = True # default pre-filter is true
    with fits.open(input_file) as f:
        data = f[1].data
        all_ra = data["RA"]
        all_dec = data["DEC"]
        all_z = data["Z"]
        colnames = data.columns.names
        all_w = data["WEIGHT"] if "WEIGHT" in colnames and use_weights else np.ones_like(all_z)
        if FKP_weights:
            all_w *= 1/(1+FKP_weights[0]*data[FKP_weights[1]]) if FKP_weights != True else data["WEIGHT_FKP"]
        if "WEIGHT" not in colnames and not FKP_weights: warn("No weights found, assigned unit weight to each particle.")
        if mask: filt = (data["STATUS"] & mask == mask) # all 1-bits from mask have to be set in STATUS; skip if mask=0
    return np.array((all_ra, all_dec, all_z, all_w)).T[filt]

def read_xi_file(xi_file: str):
    # Interpret RascalC text format using numpy functions
    if not os.path.isfile(xi_file): raise FileNotFoundError('Could not find input file %s' % xi_file)
    r_vals = np.genfromtxt(xi_file, max_rows=1)
    mu_vals = np.genfromtxt(xi_file, max_rows=1, skip_header=1)
    xi_vals = np.genfromtxt(xi_file, skip_header=2)
    return r_vals, mu_vals, xi_vals

def write_xi_file(xi_file: str, r_vals: np.ndarray[float], mu_vals: np.ndarray[float], xi_vals: np.ndarray[float]):
    # Reproduce RascalC text format using numpy functions
    header = my_a2s(r_vals) + '\n' + my_a2s(mu_vals)
    np.savetxt(xi_file, xi_vals, header=header, comments='')

def write_binning_file(out_file: str, r_edges: np.ndarray[float], print_function = blank_function):
    # Save bin edges array into a Corrfunc (and RascalC) radial binning file format
    np.savetxt(out_file, np.array((r_edges[:-1], r_edges[1:])).T)
    print_function("Binning file '%s' written successfully." % out_file)