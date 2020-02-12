### Script to convert a measured 2PCF in xi(r,mu) format to Legendre multipoles, i.e. xi_ell(r).
### This computes all even multipoles up to a specified maximum ell, approximating the integral by a sum.
### The output form is a text file with the first column specifying the r-bin, the second giving xi_0(r), the third with xi_2(r) etc.

from scipy.special import legendre
import os,sys,numpy as np

## PARAMETERS
if len(sys.argv)!=4:
    print("Usage: python convert_xi_to_multipoles.py {INFILE} {MAX_L} {OUTFILE}")
    sys.exit()

infile = str(sys.argv[1])
max_l = int(sys.argv[2])
outfile = str(sys.argv[3])

assert max_l>=0, "Maxmum multipole must be positive"
assert max_l%2==0, "Only even Legendre multipoles can be computed"
assert max_l<8, "High order multipoles cannot be reliably computed"

if not os.path.exists(infile):
    raise Exception('Could not find input file %s'%infile)

r_bins = np.genfromtxt(infile,max_rows=1)
mu_bins = np.genfromtxt(infile,max_rows=1,skip_header=1)
xi_vals = np.genfromtxt(infile,skip_header=2)

## Now convert to Legendre multipoles
xi_mult = np.zeros((len(r_bins),max_l//2+1))

for ell in np.arange(0,max_l+1,2):
    leg_mu = legendre(ell)(mu_bins) # compute mu bins

    # Compute integral as Int_0^1 dmu L_ell(mu) xi(r, mu) * (2 ell + 1)
    xi_mult[:,ell//2] = np.sum(leg_mu*xi_vals*(mu_bins[1]-mu_bins[0]),axis=1)*(2.*ell+1)

with open(outfile,"w+") as out:

    # First row contains labels
    out.write("# r-bin (Mpc/h)\t")

    for ell in np.arange(0,max_l+1,2):
        out.write("# ell = %s\t"%ell)

    # Now write data to file with each radial bin in a separate row
    for r_i,r in enumerate(r_bins):
        out.write("%.8e\t"%r)

        for ell in np.arange(0,max_l+1,2):
            out.write("%.8e\t"%xi_mult[r_i,ell//2])
        out.write("\n")

print("Output file saved to %s"%outfile)
