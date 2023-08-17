"Simple script to produce the file of mu bin Legendre factors for the C++ code in LEGENDRE_MIX mode."

from scipy.special import legendre
import numpy as np
import os, sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python mu_bin_legendre_factors.py {N_MU_BINS} {MAX_L} {OUTPUT_DIR}.")
    sys.exit(1)
n_mu_bins = int(sys.argv[1])
max_l = int(sys.argv[2])
assert max_l % 2 == 0, "Odd multipoles not supported"
n_l = max_l // 2 + 1 # number of multipoles
output_dir = str(sys.argv[3])

mu_edges = np.linspace(0, 1, n_mu_bins+1) # linearly spaced mu edges, one more than the bins

ells = np.arange(0, max_l+1, 2) # Legendre multipole indices, even only
leg_mu_factors = np.zeros((n_l, n_mu_bins))

for i, ell in enumerate(ells):
    leg_pol = legendre(ell) # Legendre polynomial
    leg_pol_int = np.polyint(leg_pol) # its indefinite integral (analytic)
    leg_mu_factors[i] = (2*ell+1)*np.diff(leg_pol_int(mu_edges)) # what we need is differences of indefinite integral between edges of mu bins = integral of Legendre polynomial over each mu bin, multiplied by 2l+1

output_file = os.path.join(output_dir, "mu_bin_legendre_factors_m%d_l%d.txt" % (n_mu_bins, n_l))
np.savetxt(output_file, leg_mu_factors.T) # more convenient to transpose so that rows are mu bins and columns are multipoles
