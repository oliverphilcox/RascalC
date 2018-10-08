## Script to generate RR pair counts from a given set of random particle. This is based on the Corrfunc code of Sinha & Garrison.
## Weights and weighted pair counts are saved in the ../weight_files/ subdirectory.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)<6:
    print "Usage: python jackknife_weights.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS}"
    sys.exit()
fname = str(sys.argv[1])
binfile = str(sys.argv[2])
mu_max = float(sys.argv[3])
nmu_bins = int(sys.argv[4])
nthreads = int(sys.argv[5])


## First read in weights and positions:
dtype = np.double 
X, Y, Z, W, J = np.genfromtxt(fname, dtype=dtype, unpack=True, usecols=[0,1,2,3,4])
J = np.array(J,dtype=int) # convert jackknives to integers
N = len(X) # number of particles
J_regions = np.unique(J) # jackknife regions in use
N_jack = len(J_regions) # number of non-empty jackknife regions

## Determine number of radial bins in binning file:
with open(binfile) as f:
    for i, l in enumerate(f):
        pass
nrbins = i + 1
print('%s radial bins are used in this file.' %nrbins)


# Now compute RR counts
from Corrfunc.theory.DDsmu import DDsmu
RR_aA=np.zeros([N_jack,nrbins*nmu_bins]);

# Iterate over jackknife regions
for i,j in enumerate(J_regions):
    filt=np.where(J==j)
    
    print("Computing pair counts for non-empty jackknife %s of %s." %(i+1,N_jack))
    
    # Compute pair counts between jackknife region and entire survey volume
    cross_RR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X,Y,Z,weights1=W,weight_type='pair_product',
                   X2=X[filt],Y2=Y[filt],Z2=Z[filt],weights2=W[filt],periodic=False,verbose=False)
    # Weight by average particle weighting
    RR_aA[i]=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
    
# Now compute weights from pair counts
w_aA=np.zeros_like(RR_aA)
RR_a = np.sum(RR_aA,axis=0)
for i,j in enumerate(J_regions):
    w_aA[i]=RR_aA[i]/RR_a # jackknife weighting for bin and region

# Save output files:
filepath = '../weight_files/'
import os
if not os.path.exists(filepath):
    os.makedirs(filepath)

weight_file='jackknife_weights_n%d_m%d_j%d.dat'%(nrbins,nmu_bins,N_jack)
with open(filepath+weight_file,"w+") as weight_file:
    for j_id,jackknife_weight in enumerate(w_aA):
        weight_file.write("%s\t" %J_regions[j_id])
        for i in range(len(jackknife_weight)):
            weight_file.write("%s\t" %jackknife_weight[i])
        weight_file.write("\n")
RR_a_file = 'binned_pair_counts_n%d_m%d_j%d.dat'%(nrbins,nmu_bins,N_jack)
with open(filepath+RR_a_file,"w+") as RR_file:
    for i in range(len(RR_a)):
        RR_file.write("%s\n" %RR_a[i])
        
print("Jackknife weights and binned pair counts written successfully to the 'weight_files/' directory")
        
