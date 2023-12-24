## Convenience function to create RR pair counts for a set of random particles, based on Corrfunc.
## This just computes the global RR pair counts, not the jackknife counts.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.

import sys
import os
import numpy as np
import math

# PARAMETERS
if len(sys.argv)!=9:
    print("Usage: python RR_counts.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} {NORMED}")
    sys.exit(1)
fname = str(sys.argv[1])
binfile = str(sys.argv[2])
mu_max = float(sys.argv[3])
nmu_bins = int(sys.argv[4])
nthreads = int(sys.argv[5])
periodic = int(sys.argv[6])
outdir=str(sys.argv[7])
normed=int(sys.argv[8])

## First read in weights and positions:

print("Reading in data");
X, Y, Z, W = np.loadtxt(fname, usecols=range(4)).T
    
N = len(X) # number of particles
weight_sum = np.sum(W)#  normalization by summed weights

print("Number of random particles %.1e"%N)

## Determine number of radial bins in binning file:
print("Counting lines in binfile");
with open(binfile) as f:
    for i, l in enumerate(f):
        pass
nrbins = i + 1
print('%s radial bins are used in this file.' %nrbins)

if not periodic:
    # Compute RR counts for the non-periodic case (measuring mu from the radial direction)
    print("Using non-periodic input data");
    def coord_transform(x,y,z):
        # Convert the X,Y,Z coordinates into Ra,Dec,comoving_distance (for use in corrfunc)
        # Shamelessly stolen from astropy
        xsq = x ** 2.
        ysq = y ** 2.
        zsq = z ** 2.

        com_dist = (xsq + ysq + zsq) ** 0.5
        s = (xsq + ysq) ** 0.5 

        if np.isscalar(x) and np.isscalar(y) and np.isscalar(z):
            Ra = math.atan2(y, x)*180./np.pi
            Dec = math.atan2(z, s)*180./np.pi
        else:
            Ra = np.arctan2(y, x)*180./np.pi+180.
            Dec = np.arctan2(z, s)*180./np.pi

        return com_dist, Ra, Dec

    # Convert coordinates to spherical coordinates
    com_dist,Ra,Dec = coord_transform(X,Y,Z);

    # Now compute RR counts
    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
    
    print("Computing pair counts")
    RR=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
                   verbose=False,is_comoving_dist=True)
    # Weight by average particle weighting
    RR_counts=RR[:]['npairs']*RR[:]['weightavg']
    if normed:
        RR_counts/=np.sum(W)**2.
        
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu
    
    print("Computing pair counts")
    RR=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,X,Y,Z,weights1=W,weight_type='pair_product',
             periodic=False,verbose=False)
    # Weight by average particle weighting
    RR_counts=RR[:]['npairs']*RR[:]['weightavg']
    if normed:
        RR_counts/=np.sum(W)**2.
    
    
# Make sure output dir exists 
if len(outdir)>0:
    os.makedirs(outdir, exist_ok=1)

outfile = os.path.join(outdir, "RR_counts_n%d_m%d_11.dat"%(nrbins,nmu_bins))
print("Saving binned pair counts as %s" %outfile);
with open(outfile,"w+") as RRfile:
    for i in range(len(RR_counts)):
        RRfile.write("%.8e\n" %RR_counts[i])
