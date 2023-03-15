### Script to smoothen a measured 2PCF in xi(r,mu) format.
### This computes all even multipoles xi_ell(r) up to a specified maximum ell, approximating the integral by a sum.
### Then r^2*xi_ell(r) are smoothened with Savitzky-Golay filter
### Finally, they are transformed back to xi(r,mu) format.

from scipy.special import legendre
from scipy.signal import savgol_filter
import os, sys, numpy as np

## PARAMETERS
if len(sys.argv)!=6:
    print("Usage: python smoothen_xi.py {INFILE} {MAX_L} {RADIAL_WINDOW_LENGTH} {RADIAL_POLYORDER} {OUTFILE}")
    sys.exit(1)

infile = str(sys.argv[1])
max_l = int(sys.argv[2])
window_length = int(sys.argv[3])
polyorder = int(sys.argv[4])
outfile = str(sys.argv[5])

assert max_l>=0, "Maxmum multipole must be positive"
assert max_l%2==0, "Only even Legendre multipoles can be computed"
assert max_l<8, "High order multipoles cannot be reliably computed"

if not os.path.exists(infile):
    raise Exception('Could not find input file %s'%infile)

r_bins = np.loadtxt(infile, max_rows=1)
mu_bins = np.loadtxt(infile, max_rows=1, skiprows=1)
xi_vals = np.loadtxt(infile, skiprows=2)

mu_edges = np.linspace(0, 1, len(mu_bins)+1) # edges of the mu bins, assumes uniform

## Now convert to Legendre multipoles
xi_mult = np.zeros((max_l//2+1, len(r_bins)))

for ell in np.arange(0, max_l+1, 2):
    leg_pol = legendre(ell) # Legendre polynomial
    leg_pol_int = np.polyint(leg_pol) # its indefinite integral (analytic)
    leg_mu_ints = np.diff(leg_pol_int(mu_edges)) # differences of indefinite integral between edges of mu bins = integral of Legendre polynomial over each mu bin

    # Compute integral as Int_0^1 dmu L_ell(mu) xi(r, mu) * (2 ell + 1)
    xi_mult[ell//2] = np.sum(leg_mu_ints * xi_vals, axis=1) * (2*ell+1)

## Perform radial smoothing
xi_mult *= r_bins**2 # multiply by r^2
xi_mult = savgol_filter(xi_mult, window_length, polyorder, axis=-1) # apply filter
xi_mult /= r_bins**2 # divide by r^2 to restore original form

## Convert back to xi(r, mu)
xi_vals = 0

for ell in np.arange(0, max_l+1, 2):
    leg_pol = legendre(ell) # Legendre polynomial
    leg_pol_int = np.polyint(leg_pol) # its indefinite integral (analytic)
    leg_mu_ints = np.diff(leg_pol_int(mu_edges)) # differences of indefinite integral between edges of mu bins = integral of Legendre polynomial over each mu bin
    leg_mu_avg = leg_mu_ints / np.diff(mu_edges) # divide integrals by bin widths to get bin-averaged value of Legendre polynomial

    xi_vals += xi_mult[ell//2][:, None] * leg_mu_avg[None, :] # accumulate xi(r, mu)

## Custom array to string function
def my_a2s(a, fmt='%.18e'):
    return ' '.join([fmt % e for e in a])

## Write to file using numpy funs
header = my_a2s(r_bins)+'\n'+my_a2s(mu_bins)
np.savetxt(outfile, xi_vals, header=header, comments='')

print("Output file saved to %s"%outfile)
