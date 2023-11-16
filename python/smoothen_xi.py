### Script to smoothen a measured 2PCF in xi(r,mu) format.
### This computes all even multipoles xi_ell(r) up to a specified maximum ell, approximating the integral by a sum.
### Then r^2*xi_ell(r) are smoothened with Savitzky-Golay filter
### Finally, they are transformed back to xi(r,mu) format.

from scipy.signal import savgol_filter
import sys, numpy as np
from mu_bin_legendre_factors import compute_mu_bin_legendre_factors
from utils import read_xi_file, write_xi_file


def smoothen_xi_multipoles(xi_mult: np.ndarray[float], r_vals: np.ndarray[float], window_length: int, polyorder: int):
    return savgol_filter(xi_mult * r_vals**2, window_length, polyorder, axis=0) / r_vals**2 # apply filter to r^2 xi_ell and divide to get back xi_ell

def smoothen_xi(xi_vals: np.ndarray[float], r_vals: np.ndarray[float], mu_edges: np.ndarray[float], max_l: int, window_length: int, polyorder: int):
    # Smoothen r, mu binned correlation function

    # Get projection factors to multipoles and back
    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse = True)
    # Project to Legendre multipoles
    xi_mult = xi_vals.dot(mu_leg_factors)
    # Perform radial smoothing
    xi_mult = smoothen_xi_multipoles(xi_mult, r_vals, window_length, polyorder)
    # Project back to r, mu binned
    xi_vals = xi_mult.dot(leg_mu_factors)

    return xi_vals

def smoothen_xi_files(infile: str, max_l: int, window_length: int, polyorder: int, outfile: str):
    r_vals, mu_vals, xi_vals = read_xi_file(infile)

    mu_edges = np.linspace(0, 1, len(mu_vals)+1) # edges of the mu bins, assumes uniform
    xi_vals_smoothed = smoothen_xi(xi_vals, r_vals, mu_edges, max_l, window_length, polyorder)

    write_xi_file(outfile, r_vals, mu_vals, xi_vals_smoothed)
    print("Output file saved to %s"%outfile)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 6:
        print("Usage: python smoothen_xi.py {INFILE} {MAX_L} {RADIAL_WINDOW_LENGTH} {RADIAL_POLYORDER} {OUTFILE}")
        sys.exit(1)

    infile = str(sys.argv[1])
    max_l = int(sys.argv[2])
    window_length = int(sys.argv[3])
    polyorder = int(sys.argv[4])
    outfile = str(sys.argv[5])

    smoothen_xi_files(infile, max_l, window_length, polyorder, outfile)