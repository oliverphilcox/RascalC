## Script to generate RR pair counts from two distinct sets of random particles. This is based on the Corrfunc code of Sinha & Garrison.
## Here we compute the cross-weight RR_{XY} between two sets of particles.
## Weights and weighted pair counts are saved in the {OUTPUT_DIR} subdirectory. 
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=9:
    print("Usage: python jackknife_weights_cross.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR}")
    sys.exit()
fname = str(sys.argv[1])
fname2 = str(sys.argv[2])
binfile = str(sys.argv[3])
mu_max = float(sys.argv[4])
nmu_bins = int(sys.argv[5])
nthreads = int(sys.argv[6])
periodic = int(sys.argv[7])
outdir=str(sys.argv[8])


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
weight_sum = np.sum(W)
weight_sum2 = np.sum(W2)
J_regions = np.unique(np.concatenate([J,J2])) # jackknife regions in use
J_regions.sort() # sort to ensure same ordering everywhere
N_jack = len(J_regions) # number of jackknife regions which are non-empty in at least one sample

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

    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
    
    def non_periodic_RR_counts(Ra,Dec,com_dist,W,Ra2,Dec2,com_dist2,W2,distinct=False):
        # General function to compute RR pair counts for two possibly distinct fields
        RR_aA=np.zeros([N_jack,nrbins*nmu_bins]);
        
        # Iterate over jackknife regions
        for i,j in enumerate(J_regions):
            # Compute pair counts between jackknife region and entire survey volume
            print("Computing pair counts for non-empty jackknife %s of %s." %(i+1,N_jack))
            
            # First compute contribution from JK region in data1 and all of data2:
            filt=np.where(J==j)
            if len(filt[0])>0:
                cross_RR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,Ra2,Dec2,com_dist2,weights1=W2,weight_type='pair_product',
                            RA2=Ra[filt],DEC2=Dec[filt],CZ2=com_dist[filt],weights2=W[filt],verbose=False,is_comoving_dist=True)
                # Weight by average particle weighting
                RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
            
            if distinct:
                # This term is only needed if field1!=field2 else it arises from double counting
                # Now compute contribution from JK region in data2 and all of data1:
                filt2 = np.where(J2==j)
                if len(filt2[0])>0:
                    cross_RR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
                                RA2=Ra2[filt2],DEC2=Dec2[filt2],CZ2=com_dist2[filt2],weights2=W2[filt2],verbose=False,is_comoving_dist=True)
                    # Weight by average particle weighting
                    RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                RR_aA[i]/=2.
        return RR_aA
    
    # Now compute RR pair counts
    print("Computing Field 1 x Field 1 pair counts, w_aA^{11}")
    RR_aA_11=non_periodic_RR_counts(Ra,Dec,com_dist,W,Ra,Dec,com_dist,W,distinct=False)
    print("Computing Field 1 x Field 2 pair counts, w_aA^{12}")
    RR_aA_12=non_periodic_RR_counts(Ra,Dec,com_dist,W,Ra2,Dec2,com_dist2,W2,distinct=True)
    print("Computing Field 2 x Field 2 pair counts, w_aA^{22}")
    RR_aA_22=non_periodic_RR_counts(Ra2,Dec2,com_dist2,W2,Ra2,Dec2,com_dist2,W2,distinct=False)
    
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu

    def periodic_RR_counts(X,Y,Z,W,X2,Y2,Z2,W2,distinct=False):
        # Compute the RR pair counts for two possibly distinct fields
        RR_aA=np.zeros([N_jack,nrbins*nmu_bins]);

        # Iterate over jackknife regions
        for i,j in enumerate(J_regions):
            print("Computing pair counts for non-empty jackknife %s of %s." %(i+1,N_jack))
            # Compute pair counts between jackknife region and entire survey volume
            
            # First compute contribution from JK region in data1 and all of data2:
            filt=np.where(J==j)
            print(filt,len(filt))
            if len(filt[0])>0:
                cross_RR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X2,Y2,Z2,weights1=W2,weight_type='pair_product',
                            X2=X[filt],Y2=Y[filt],Z2=Z[filt],weights2=W[filt],periodic=True,verbose=False)
                # Weight by average particle weighting
                RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
            
            if distinct:
                # Now compute contribution from JK region in data2 and all of data1:
                filt2 = np.where(J2==j)
                if len(filt2[0])>0:
                    cross_RR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,X,Y,Z,weights1=W,weight_type='pair_product',
                                X2=X2[filt],Y2=Y2[filt],Z2=Z2[filt],weights2=W2[filt],periodic=True,verbose=False)
                    # Weight by average particle weighting
                    RR_aA[i]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                RR_aA[i]/=2.
            return RR_aA

    # Now compute RR pair counts
    print("Computing Field 1 x Field 1 pair counts, w_aA^{11}")
    RR_aA_11=periodic_RR_counts(X,Y,Z,W,X,Y,Z,W,distinct=False)
    print("Computing Field 1 x Field 2 pair counts, w_aA^{12}")
    RR_aA_12=periodic_RR_counts(X,Y,Z,W,X2,Y2,Z2,W2,distinct=True)
    print("Computing Field 2 x Field 2 pair counts, w_aA^{22}")
    RR_aA_22=periodic_RR_counts(X2,Y2,Z2,W2,X2,Y2,Z2,W2,distinct=False)
  
def summed_weights(RR_aA):
    # Now compute weights from pair counts
    w_aA=np.zeros_like(RR_aA)
    RR_a = np.sum(RR_aA,axis=0)
    for i,j in enumerate(J_regions):
        w_aA[i]=RR_aA[i]/RR_a # jackknife weighting for bin and region
    return w_aA,RR_a

w_aA_11,RR_a_11 = summed_weights(RR_aA_11)
w_aA_12,RR_a_12 = summed_weights(RR_aA_12)
w_aA_22,RR_a_22 = summed_weights(RR_aA_22)

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

all_weights = [w_aA_11,w_aA_12,w_aA_22]
all_pairs = [RR_a_11,RR_a_12,RR_a_22]
all_counts = [RR_aA_11,RR_aA_12,RR_aA_22]
all_indices = ['11','12','22']

for index in range(3):
    weight_file='jackknife_weights_n%d_m%d_j%d_%s.dat'%(nrbins,nmu_bins,N_jack,all_indices[index])
    print("Saving jackknife weight as %s"%weight_file)
    with open(os.path.join(outdir, weight_file), "w+") as weight_file:
        for j_id,jackknife_weight in enumerate(all_weights[index]):
            weight_file.write("%s\t" %J_regions[j_id])
            for i in range(len(jackknife_weight)):
                weight_file.write("%s" %jackknife_weight[i])
                if i == len(jackknife_weight)-1:
                    weight_file.write("\n");
                else:
                    weight_file.write("\t");
    RR_a_file = 'binned_pair_counts_n%d_m%d_j%d_%s.dat'%(nrbins,nmu_bins,N_jack,all_indices[index])
    print("Saving binned pair counts as %s" %RR_a_file);
    with open(os.path.join(outdir, RR_a_file), "w+") as RR_file:
        for i in range(len(all_pairs[index])):
            RR_file.write("%s\n" %all_pairs[index][i])
    
    RR_aA_file = 'jackknife_pair_counts_n%d_m%d_j%d_%s.dat'%(nrbins,nmu_bins,N_jack,all_indices[index])
    print("Saving normalized jackknife pair counts as %s"%RR_aA_file)
    with open(os.path.join(outdir, RR_aA_file), "w+") as jackRR_file:
        for j_id,pair_count in enumerate(all_counts[index]):
            this_jk = J_regions[j_id]
            if index==0:
                norm = weight_sum**2.#np.sum(W[J==this_jk)
            if index==1:
                norm = weight_sum*weight_sum2#np.sum(W2[J2==this_jk])+weight_sum2*np.sum(W[J==this_jk])
            if index==2:
                norm = weight_sum2**2.#*np.sum(W2[J2==this_jk])
            jackRR_file.write("%d\t" %this_jk)
            for i in range(len(pair_count)):
                jackRR_file.write("%.8e" %(pair_count[i]/norm))
                if i == len(pair_count)-1:
                    jackRR_file.write("\n");
                else:
                    jackRR_file.write("\t");
        
print("Jackknife weights and binned pair counts written successfully to the %s directory"%outdir)
        
