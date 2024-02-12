from .utils import blank_function
import numpy as np


def write_binning_file(out_file: str, r_edges: np.ndarray[float], print_function = blank_function):
    # Save bin edges array into a Corrfunc (and RascalC) radial binning file format
    np.savetxt(out_file, np.array((r_edges[:-1], r_edges[1:])).T)
    print_function("Binning file '%s' written successfully." % out_file)


def write_binning_file_linear(out_file: str, r_min: float, r_max: float, nrbins: int, print_function = print):
    print_function("Using LINEAR binning")
    # Define radial bins
    r_edges = np.linspace(r_min, r_max, nrbins+1)
    if r_edges[0] <= 1e-4: r_edges[0] = 1e-4 # exclude very small separations for stability
    write_binning_file(out_file, r_edges, print_function)


def write_binning_file_log(out_file: str, r_min: float, r_max: float, nrbins: int, print_function = print):
    print_function("Using LOG binning")
    if r_min <= 0: raise ValueError('Minimum r must be positive to take logarithm')
    # Define radial bins
    r_edges = np.geomspace(r_min, r_max, nrbins+1)
    write_binning_file(out_file, r_edges, print_function)


def write_binning_file_hybrid(out_file: str, r_min: float, r_cut: float, r_max: float, n_log_bins: int, n_lin_bins: int, print_function = print):
    print_function("Using hybrid binning: log binning up to r = %.1f, linear binning up to r = %.1f" %(r_cut,r_max))
    if r_min <= 0: raise ValueError('Minimum r must be positive to take logarithm')
    # Define radial bins
    r_edges = np.concatenate((np.geomspace(r_min, r_cut, n_log_bins+1)[:-1], np.linspace(r_cut, r_max, n_lin_bins+1)))
    write_binning_file(out_file, r_edges, print_function)