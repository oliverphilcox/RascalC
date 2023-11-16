# Contains some utility functions widely used in other scripts
# Not intended for execution from command line
import sys, os
import numpy as np

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

def symmetrized(A):
    # symmetrize a 2D matrix
    return 0.5 * (A + A.T)

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