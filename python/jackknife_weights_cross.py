## Script to generate RR pair counts from two sets of random particles. This is based on the Corrfunc code of Sinha & Garrison.
## Here we compute the cross-weight RR_{XY} between two sets of particles.
## Weights and weighted pair counts are saved in the {OUTPUT_DIR} subdirectory. 
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.


import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=11:
    print("Usage: python jackknife_weights.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {MU_MAX} {N_MU_BINS} {N_GAL_1} {N_GAL_2} {NTHREADS} {PERIODIC} {OUTPUT_DIR}")
    sys.exit()
fname = str(sys.argv[1])
fname2 = str(sys.argv[2])
binfile = str(sys.argv[3])
mu_max = float(sys.argv[4])
nmu_bins = int(sys.argv[5])
N_gal = float(sys.argv[6])
N_gal2 = float(sys.argv[7])
nthreads = int(sys.argv[8])
periodic = int(sys.argv[9])
outdir=str(sys.argv[10])


## First read in weights and positions:
dtype = np.double 

print("Counting lines in file 1 of 2")
total_lines=0
for n, line in enumerate(open(fname, 'r')):
    total_lines+=1
print("Counting lines in file 2 of 2")
total_lines2=0
for n, line in enumerate(open(fname2,'r')):
    total_lines2+=1

X,Y,Z,W,J=[np.zeros(total_lines) for _ in range(5)]
X2,Y2,Z2,W2,J2=[np.zeros(total_lines2) for _ in range(5)]

print("\nReading in data from file 1:");
for n, line in enumerate(open(fname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    X[n]=split_line[0];
    Y[n]=split_line[1];
    Z[n]=split_line[2];
    W[n]=split_line[3];
    J[n]=int(split_line[4]);

print("\nReading in data from file 2:");
for n, line in enumerate(open(fname2, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines2))
    split_line=np.array(line.split(" "), dtype=float) 
    X2[n]=split_line[0];
    Y2[n]=split_line[1];
    Z2[n]=split_line[2];
    W2[n]=split_line[3];
    J2[n]=int(split_line[4]);

N = len(X) # number of particles
N2 = len(X2)
J_regions = np.unique(np.concatenate([J,J2])) # jackknife regions in use
N_jack = len(J_regions) # number of jackknife regions which are non-empty in at least one sample

print("Ratio of number of galaxies to number of random particles: %.1e (dataset 1) %.1e (dataset 2)" %(N_gal/N,N_gal2/N2))

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
    com_dist2,Ra2,Dec2 = coord_transform(X2,Y2,Z2);

    # Now compute RR counts
    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
    RR_aA=np.zeros([N_jack,nrbins*nmu_bins]);

    # Iterate over jackknife regions
    for i,j in enumerate(J_regions):
        # Compute pair counts between jackknife region and entire survey volume
        print("Computing pair counts for non-empty jackknife %s of %s." %(i+1,N_jack))
        
        # First compute contribution from JK region in data1 and all of data2:
        filt=np.where(J==j)
        if len(filt)>0:
            cross_RR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,Ra2,Dec2,com_dist2,weights1=W2,weight_type='pair_product',
                        RA2=Ra[filt],DEC2=Dec[filt],CZ2=com_dist[filt],weights2=W[filt],verbose=False,is_comoving_dist=True)
            # Weight by average particle weighting
            RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']*(N_gal/N)*(N_gal2/N2)
        
        # Now compute contribution from JK region in data2 and all of data1:
        filt2 = np.where(J2==j)
        if len(filt2)>0:
            cross_RR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
                        RA2=Ra2[filt2],DEC2=Dec2[filt2],CZ2=com_dist2[filt2],weights2=W2[filt2],verbose=False,is_comoving_dist=True)
            # Weight by average particle weighting
            RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']*(N_gal/N)*(N_gal2/N2)
        
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu
    
    RR_aA=np.zeros([N_jack,nrbins*nmu_bins]);

    # Iterate over jackknife regions
    for i,j in enumerate(J_regions):
        print("Computing pair counts for non-empty jackknife %s of %s." %(i+1,N_jack))
        # Compute pair counts between jackknife region and entire survey volume
        
        # First compute contribution from JK region in data1 and all of data2:
        filt=np.where(J==j)
        if len(filt)>0:
            cross_RR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X2,Y2,Z2,weights1=W2,weight_type='pair_product',
                           X2=X[filt],Y2=Y[filt],Z2=Z[filt],weights2=W[filt],periodic=False,verbose=False)
            # Weight by average particle weighting
            RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']*(N_gal/N)*(N_gal2/N2)
        
        # Now compute contribution from JK region in data2 and all of data1:
        filt2 = np.where(J2==j)
        if len(filt2)>0:
            cross_RR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X,Y,Z,weights1=W,weight_type='pair_product',
                           X2=X2[filt],Y2=Y2[filt],Z2=Z2[filt],weights2=W2[filt],periodic=False,verbose=False)
            # Weight by average particle weighting
            RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']*(N_gal/N)*(N_gal2/N2)

# Now compute weights from pair counts
w_aA=np.zeros_like(RR_aA)
RR_a = np.sum(RR_aA,axis=0)
for i,j in enumerate(J_regions):
    w_aA[i]=RR_aA[i]/RR_a # jackknife weighting for bin and region

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

weight_file='jackknife_weights_n%d_m%d_j%d.dat'%(nrbins,nmu_bins,N_jack)
print("Saving jackknife weight as %s"%weight_file)
with open(outdir+weight_file,"w+") as weight_file:
    for j_id,jackknife_weight in enumerate(w_aA):
        weight_file.write("%s\t" %J_regions[j_id])
        for i in range(len(jackknife_weight)):
            weight_file.write("%s" %jackknife_weight[i])
            if i == len(jackknife_weight)-1:
                weight_file.write("\n");
            else:
                weight_file.write("\t");
RR_a_file = 'binned_pair_counts_n%d_m%d_j%d.dat'%(nrbins,nmu_bins,N_jack)
print("Saving binned pair counts as %s" %RR_a_file);
with open(outdir+RR_a_file,"w+") as RR_file:
    for i in range(len(RR_a)):
        RR_file.write("%s\n" %RR_a[i])
        
print("Jackknife weights and binned pair counts written successfully to the %s directory"%outdir)
        
