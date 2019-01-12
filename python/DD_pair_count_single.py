## Script to generate DD pair counts from a given set of galaxies particles. This is based on the Corrfunc code of Sinha & Garrison.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=8:
    print("Usage: python DD_pair_counts.py {DATA_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR}")
    sys.exit()
fname = str(sys.argv[1])
binfile = str(sys.argv[2])
mu_max = float(sys.argv[3])
nmu_bins = int(sys.argv[4])
nthreads = int(sys.argv[5])
periodic = int(sys.argv[6])
outdir=str(sys.argv[7])


## First read in weights and positions:
dtype = np.double 

print("Counting lines in file")
total_lines=0
for n, line in enumerate(open(fname, 'r')):
    total_lines+=1

X,Y,Z,W=[np.zeros(total_lines) for _ in range(4)]

print("Reading in data");
for n, line in enumerate(open(fname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    X[n]=split_line[0];
    Y[n]=split_line[1];
    Z[n]=split_line[2];
    W[n]=split_line[3];
    
N = len(X) # number of particles

print("Number of galaxy particles %.1e"%N)

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
    
    all_DD=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',verbose=True,is_comoving_dist=True)
    # Weight by average particle weighting
    DD_aA=all_DD[:]['npairs']*all_DD[:]['weightavg']#*(N_gal/N)**2.
        
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu
    
    all_DD=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X,Y,Z,weights1=W,weight_type='pair_product',periodic=False,verbose=True)
    # Weight by average particle weighting
    DD_aA[i]=all_DD[:]['npairs']*all_DD[:]['weightavg']
    

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

outfile='DD_n%d_m%d.dat'%(nrbins,nmu_bins)
print("Saving DD pair counts as %s"%outfile)
with open(outdir+outfile,"w+") as outfile:
    for i in range(len(DD_aA)):
        outfile.write("%s\n" %DD_aA[i])
        
print("DD pair counts written successfully to the %s directory"%outdir)
        
