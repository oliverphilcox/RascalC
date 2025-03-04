import numpy as np
import os
from ..utils import my_a2s


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