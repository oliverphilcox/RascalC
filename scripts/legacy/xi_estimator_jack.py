## Script to compute jackknife estimates of the correlation function xi^J(r,mu) via the Landay-Szelay estimator for a single set of random particles.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.
## This must be binned in the same binning as the desired covariance matrix

import sys
import numpy as np
import math

# PARAMETERS
if len(sys.argv)!=10:
    if len(sys.argv)!=11:
        print("Usage: python xi_estimator_jack.py {GALAXY_FILE} {RANDOM_FILE_DR} {RANDOM_FILE_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_jackknife_counts}]")
        sys.exit(1)
Dname = str(sys.argv[1])
RnameDR = str(sys.argv[2])
RnameRR = str(sys.argv[3])
binfile = str(sys.argv[4])
mu_max = float(sys.argv[5])
nmu_bins = int(sys.argv[6])
nthreads = int(sys.argv[7])
periodic = int(sys.argv[8])
outdir=str(sys.argv[9])

if len(sys.argv)==11:
    print("Using pre-defined RR counts")
    RRname=str(sys.argv[10])
else:
    RRname=""

## First read in weights and positions:
# Read first set of randoms

print("Reading in DR random data")
rX_DR, rY_DR, rZ_DR, rW_DR, rJ_DR = np.loadtxt(RnameDR, usecols=range(5)).T
rJ_DR = rJ_DR.astype(int) # jackknife region is integer

N_randDR = len(rX_DR) # number of particles

if len(RRname)==0:
    if RnameRR!=RnameDR:
        # only read in RR file if distinct from DR and needed by the code:

        print("Reading in RR random data")
        rX_RR, rY_RR, rZ_RR, rW_RR, rJ_RR = np.loadtxt(RnameRR, usecols=range(5)).T
        rJ_RR = rJ_RR.astype(int) # jackknife region is integer
            
        N_randRR = len(rX_RR) # number of particles
    else:
        # just copy if its the same
        rX_RR=rX_DR
        rY_RR=rY_DR
        rZ_RR=rZ_DR
        rW_RR=rW_DR
        rJ_RR=rJ_DR
        N_randRR=N_randDR
else:
    # empty placeholders
    rX_RR,rY_RR,rZ_RR,rW_RR,rJ_RR=[],[],[],[],[]
    print("Counting lines in RR random file")
    N_randRR=0
    for n, line in enumerate(open(RnameRR, 'r')):
        N_randRR+=1

print("Reading in galaxy data")
dX, dY, dZ, dW, dJ = np.loadtxt(Dname, usecols=range(5)).T
dJ = dJ.astype(int) # jackknife region is integer

N_gal = len(dX) # number of particles

print("Number of random particles %.1e (DR) %.1e (RR)"%(N_randDR, N_randRR))
print("Number of galaxy particles %.1e"%N_gal)

# Determine number of jackknives
if len(RRname)==0:
    J_regions = np.unique(np.concatenate([rJ_RR,rJ_DR,dJ]))
else:
    J_regions = np.unique(np.concatenate([rJ_DR,dJ]))

N_jack = len(J_regions)

print("Using %d non-empty jackknife regions"%N_jack)

## Determine number of radial bins in binning file:
print("Counting lines in binfile");
with open(binfile) as f:
    for i, l in enumerate(f):
        pass
nrbins = i + 1
all_bins = np.loadtxt(binfile)
mean_bins=0.5*(all_bins[:,0]+all_bins[:,1])
print('%s radial bins are used in this file in the range [%d, %d]' %(nrbins,all_bins[0,0],all_bins[-1,1]))

if not periodic:
    # Compute RR, DR and DD counts for the non-periodic case (measuring mu from the radial direction)
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
    r_com_dist_DR,r_Ra_DR,r_Dec_DR = coord_transform(rX_DR,rY_DR,rZ_DR);
    d_com_dist,d_Ra,d_Dec = coord_transform(dX,dY,dZ);

    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
    
    import time
    init=time.time()
    
    # Now compute RR counts
    RR_counts = np.zeros([N_jack,nrbins*nmu_bins])    
    if len(RRname)!=0:
        RRfile = np.loadtxt(RRname) # read pre-computed RR counts
        if len(RRfile[0,:])!=(1+nrbins*nmu_bins):
            raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        if len(RR_counts[:,0])!=N_jack:
            raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        for jk in range(N_jack):
            RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
            # NB: these are already normalized
    else:
        r_com_dist_RR,r_Ra_RR,r_Dec_RR = coord_transform(rX_RR,rY_RR,rZ_RR);
        # Compute RR pair counts
        for i,j in enumerate(J_regions):
            # Compute pair counts between jackknife region and entire survey regions
            filt = np.where(rJ_RR==j)
            print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if len(filt[0])>0:
                cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,
                                      weight_type='pair_product',RA2=r_Ra_RR[filt],DEC2=r_Dec_RR[filt],CZ2=r_com_dist_RR[filt],
                                      weights2=rW_RR[filt],verbose=False,is_comoving_dist=True)
                RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                RR_counts[i,:] /= np.sum(rW_RR)**2. # normalize by product of sum of weights
    print("Finished RR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DR counts
    DR_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        
        # Compute pair counts between data jackknife and random survey
        filt = np.where(dJ==j)
        if len(filt[0])>0:
            cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],
                                   weights1=dW[filt],weight_type='pair_product', RA2=r_Ra_DR, DEC2=r_Dec_DR, 
                                   CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=False,is_comoving_dist=True)
            DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
        
        # Compute pair coutnts between random jackknife and data survey
        filt2 = np.where(rJ_DR==j)
        if len(filt2[0])>0:
            cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                   weights1=dW,weight_type='pair_product', RA2=r_Ra_DR[filt2], DEC2=r_Dec_DR[filt2], 
                                   CZ2 = r_com_dist_DR[filt2], weights2 = rW_DR[filt2], verbose=False,is_comoving_dist=True)
            DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
        DR_counts[i,:] /= np.sum(rW_DR)*np.sum(dW) # normalize by product of sum of weights
    print("Finished DR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DD counts
    DD_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        # Compute pair counts between jackknife region and entire survey regions
        filt = np.where(dJ==j)
        print("Computing DD pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        if len(filt[0])>0:
            cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                                    weight_type='pair_product',RA2=d_Ra[filt],DEC2=d_Dec[filt],CZ2=d_com_dist[filt],
                                    weights2=dW[filt],verbose=False,is_comoving_dist=True)
            DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
            DD_counts[i,:] /= np.sum(dW)**2. # normalize by product of sum of weights
    print("Finished after %d seconds"%(time.time()-init))
    
else:
    # Compute xi for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu
    
    import time
    init = time.time()
    
    # Now compute RR counts
    RR_counts = np.zeros([N_jack,nrbins*nmu_bins])    
    if len(RRname)!=0:
        RRfile = np.loadtxt(RRname) # read pre-computed RR counts
        if len(RRfile[0,:])!=(1+nrbins*nmu_bins):
            raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        if len(RR_counts[:,0])!=N_jack:
            raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        for jk in range(N_jack):
            RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
            # NB: these are already normalized 
    else:
        # Compute RR pair counts
        for i,j in enumerate(J_regions):
            # Compute pair counts between jackknife region and entire survey regions
            filt = np.where(rJ_RR==j)
            print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if len(filt[0])>0:
                cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX_RR,rY_RR,rZ_RR,weights1=rW_RR,
                                      weight_type='pair_product',X2=rX_RR[filt],Y2=rY_RR[filt],Z2=rZ_RR[filt],
                                      weights2=rW_RR[filt],verbose=False,periodic=True)
                RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                RR_counts[i,:] /= np.sum(rW_RR)**2. # normalize by product of sum of weights
    print("Finished RR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DR counts
    DR_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        
        # Compute pair counts between data jackknife and random survey
        filt = np.where(dJ==j)
        if len(filt[0])>0:
            cross_DR1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],
                                   weights1=dW[filt],weight_type='pair_product', X2=rX_DR, Y2=rY_DR, 
                                   Z2 = rZ_DR, weights2 = rW_DR, verbose=False,periodic=True)
            DR_counts[i,:] += 0.5*cross_DR1[:]['npairs']*cross_DR1[:]['weightavg']
        
        # Compute pair coutnts between random jackknife and data survey
        filt2 = np.where(rJ_DR==j)
        if len(filt2[0])>0:
            cross_DR2 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                   weights1=dW,weight_type='pair_product', X2=rX_DR[filt2], Y2=rY_DR[filt2], 
                                   Z2 = rZ_DR[filt2], weights2 = rW_DR[filt2], verbose=False,periodic=True)
            DR_counts[i,:] += 0.5*cross_DR2[:]['npairs']*cross_DR2[:]['weightavg']
        DR_counts[i,:] /= np.sum(rW_DR)*np.sum(dW) # normalize by product of sum of weights
    print("Finished DR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DD counts
    DD_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        # Compute pair counts between jackknife region and entire survey regions
        filt = np.where(dJ==j)
        print("Computing DD pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        if len(filt[0])>0:
            cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,
                                    weight_type='pair_product',X2=dX[filt],Y2=dY[filt],Z2=dZ[filt],
                                    weights2=dW[filt],verbose=False,periodic=True)
            DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
            DD_counts[i,:] /= np.sum(dW)**2. # normalize by product of sum of weights
    print("Finished after %d seconds"%(time.time()-init))
    

# Now compute correlation function
xi_function = np.zeros_like(RR_counts)
for j in range(N_jack):
    xi_function[j] = DD_counts[j]/RR_counts[j] - 2.*DR_counts[j]/RR_counts[j] + 1.

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

outname='xi_jack_n%d_m%d_j%d_11.dat'%(nrbins,nmu_bins,N_jack)
print("Saving correlation function")
with open(os.path.join(outdir, outname), "w+") as outfile:
    for r in mean_bins:
        outfile.write("%.8e "%r)
    outfile.write("\n")
    for mu in mean_mus:
        outfile.write("%.8e "%mu)
    outfile.write("\n");
    for i in range(N_jack):
        for j in range(nrbins*nmu_bins):
            outfile.write("%.8e "%xi_function[i,j])
        outfile.write("\n")
        
print("Correlation function written successfully to %s"%(outdir+outname))

print("NB: Number of galaxies is %d"%N_gal)

#print("ADD IN NORM")

np.savez("test_jack_xi.npz",xi_jack=xi_function,DD=DD_counts,RR=RR_counts,DR=DR_counts);
