"Simple script to produce the file of mu bin Legendre factors for the C++ code in LEGENDRE_MIX mode."

from scipy.special import legendre
import numpy as np
import os, sys


def compute_mu_bin_legendre_factors(mu_edges: np.ndarray[float], max_l: int, do_inverse: bool = False) -> np.ndarray[float] | (np.ndarray[float], np.ndarray[float]):
    if max_l < 0: raise ValueError("Maximum multipole must be positive")
    if max_l % 2 != 0: raise ValueError("Odd multipoles not supported")
    n_l = max_l // 2 + 1 # number of multipoles

    ells = np.arange(0, max_l+1, 2) # Legendre multipole indices, even only
    mu_leg_factors = np.zeros((len(mu_edges), n_l))
    if do_inverse: leg_mu_factors = np.zeros((n_l, len(mu_edges)))

    for i, ell in enumerate(ells):
        leg_pol = legendre(ell) # Legendre polynomial
        leg_pol_int = np.polyint(leg_pol) # its indefinite integral (analytic)
        mu_leg_factors[:, i] = (2*ell+1)*np.diff(leg_pol_int(mu_edges)) # differences of indefinite integral between edges of mu bins = integrals of Legendre polynomial over each mu bin, multiplied by 2l+1 convert from mu-binned values to Legendre multipoles in each radial bin
        if do_inverse: leg_mu_factors[i] = np.diff(leg_pol_int(mu_edges)) / np.diff(mu_edges) # the bin-averaged values of the multipole is an approximation for inverse conversion, from Legendre multipoles to mu-binned values in each radial bin
    
    if do_inverse:
        return mu_leg_factors, leg_mu_factors
    return mu_leg_factors

def write_mu_bin_legendre_factors(n_mu_bins, max_l, output_dir):
    mu_edges = np.linspace(0, 1, n_mu_bins+1) # edges of the mu bins, assumes uniform
    mu_bin_legendre_factors = compute_mu_bin_legendre_factors(mu_edges, max_l)

    output_file = os.path.join(output_dir, "mu_bin_legendre_factors_m%d_l%d.txt" % (n_mu_bins, max_l))
    os.makedirs(output_dir, exist_ok=1) # make sure the directory exists
    np.savetxt(output_file, mu_bin_legendre_factors)
    return output_file

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python mu_bin_legendre_factors.py {N_MU_BINS} {MAX_L} {OUTPUT_DIR}.")
        sys.exit(1)
    n_mu_bins = int(sys.argv[1])
    max_l = int(sys.argv[2])
    output_dir = str(sys.argv[3])

    write_mu_bin_legendre_factors(n_mu_bins, max_l, output_dir)