### Script to convert a measured 2PCF in xi(r,mu) format to Legendre multipoles, i.e. xi_ell(r).
### This computes all even multipoles up to a specified maximum ell, approximating the integral by a sum.
### The output form is a text file with the first column specifying the r-bin, the second giving xi_0(r), the third with xi_2(r) etc.

import sys
import numpy as np
from .utils import read_xi_file
from .mu_bin_legendre_factors import compute_mu_bin_legendre_factors


def convert_xi_to_multipoles(xi_vals: np.ndarray[float], mu_edges: np.ndarray[float], max_l: int):
    ## Convert s, mu-binned correlation function to Legendre multipoles
    return xi_vals.dot(compute_mu_bin_legendre_factors(mu_edges, max_l))

def convert_xi_to_multipoles_files(infile, max_l, outfile):
    r_vals, mu_vals, xi_vals = read_xi_file(infile)

    mu_edges = np.linspace(0, 1, len(mu_vals)+1) # edges of the mu bins, assumes uniform

    xi_mult = convert_xi_to_multipoles(xi_vals, mu_edges, max_l)

    with open(outfile,"w+") as out:

        # First row contains labels
        out.write("# r-bin (Mpc/h)\t")

        for ell in np.arange(0,max_l+1,2):
            out.write("# ell = %s\t"%ell)
        out.write("\n")

        # Now write data to file with each radial bin in a separate row
        for r_i,r in enumerate(r_vals):
            out.write("%.8e\t"%r)

            for ell in np.arange(0,max_l+1,2):
                out.write("%.8e\t"%xi_mult[r_i,ell//2])
            out.write("\n")

    print("Output file saved to %s"%outfile)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 4:
        print("Usage: python convert_xi_to_multipoles.py {INFILE} {MAX_L} {OUTFILE}")
        sys.exit(1)

    infile = str(sys.argv[1])
    max_l = int(sys.argv[2])
    outfile = str(sys.argv[3])

    convert_xi_to_multipoles_files(infile, max_l, outfile)