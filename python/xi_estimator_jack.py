## Script to compute jackknife estimates of the correlation function xi^J(r,mu) via the Landay-Szelay estimator for a single set of random particles.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.
## This must be binned in the same binning as the desired covariance matrix

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=9:
    if len(sys.argv)!=10:
        print("Usage: python xi_estimator.py {GALAXY_FILE} {RANDOM_FILE} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_jackknife_counts}]")
        sys.exit()
Dname = str(sys.argv[1])
Rname = str(sys.argv[2])
binfile = str(sys.argv[3])
mu_max = float(sys.argv[4])
nmu_bins = int(sys.argv[5])
nthreads = int(sys.argv[6])
periodic = int(sys.argv[7])
outdir=str(sys.argv[8])

if len(sys.argv)==10:
    print("Using pre-defined RR counts")
    RRname=str(sys.argv[9])
else:
    RRname=""

## First read in weights and positions:
dtype = np.double 

print("Counting lines in random file")
total_lines=0
for n, line in enumerate(open(Rname, 'r')):
    total_lines+=1

rX,rY,rZ,rW,rJ=[np.zeros(total_lines) for _ in range(5)]

print("Reading in random data");
for n, line in enumerate(open(Rname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    rX[n]=split_line[0];
    rY[n]=split_line[1];
    rZ[n]=split_line[2];
    rW[n]=split_line[3];
    rJ[n]=int(split_line[4]);

N_rand = len(rX) # number of particles

print("Counting lines in galaxy file")
total_lines=0
for n, line in enumerate(open(Dname, 'r')):
    total_lines+=1

dX,dY,dZ,dW,dJ=[np.zeros(total_lines) for _ in range(5)]

print("Reading in galaxy data");
for n, line in enumerate(open(Dname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    dX[n]=split_line[0];
    dY[n]=split_line[1];
    dZ[n]=split_line[2];
    dW[n]=split_line[3];
    dJ[n]=int(split_line[4]);

N_gal = len(dX) # number of particles

print("Number of random particles %.1e"%N_rand)
print("Number of galaxies particles %.1e"%N_gal)

# Determine number of jackknifes
J_regions = np.unique(np.concatenate([rJ,dJ]))
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
    r_com_dist,r_Ra,r_Dec = coord_transform(rX,rY,rZ);
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
        if len(RR_counts[:,0])!=N_jack):
            raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        for jk in range(N_jack):
            RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
    else:
        # Compute RR pair counts
        for i,j in enumerate(J_regions):
            # Compute pair counts between jackknife region and entire survey regions
            filt = np.where(rJ==j)
            print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if len(filt[0])>0:
                cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra,r_Dec,r_com_dist,weights1=rW,
                                      weight_type='pair_product',RA2=r_Ra[filt],DEC2=r_Dec[filt],CZ2=r_com_dist[filt],
                                      weights2=rW[filt],verbose=False,is_comoving_dist=True)
                RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
    print("Finished RR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DR counts
    DR_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        
        # Compute pair counts between data jackknife and random survey
        filt = np.where(dJ==j)
        if len(filt[0])>0:
            cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],
                                   weights1=dW[filt],weight_type='pair_product', RA2=r_Ra, DEC2=r_Dec, 
                                   CZ2 = r_com_dist, weights2 = rW, verbose=False,is_comoving_dist=True)
            DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
        
        # Compute pair coutnts between random jackknife and data survey
        filt = np.where(rJ==j)
        if len(filt[0])>0:
            cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                   weights1=dW,weight_type='pair_product', RA2=r_Ra[filt], DEC2=r_Dec[filt], 
                                   CZ2 = r_com_dist[filt], weights2 = rW[filt], verbose=False,is_comoving_dist=True)
            DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
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
    print("Finished after %d seconds"%(time.time()-init))
    
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
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
        if len(RR_counts[:,0])!=N_jack):
            raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
        for jk in range(N_jack):
            RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
    else:
        # Compute RR pair counts
        for i,j in enumerate(J_regions):
            # Compute pair counts between jackknife region and entire survey regions
            filt = np.where(rJ==j)
            print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if len(filt[0])>0:
                cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX,rY,rZ,weights1=rW,
                                      weight_type='pair_product',X2=rX[filt],Y2=rY[filt],Z2=rZ[filt],
                                      weights2=rW[filt],verbose=False,periodic=True)
                RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
    print("Finished RR pair counts after %d seconds"%(time.time()-init))
    
    # Now compute DR counts
    DR_counts = np.zeros_like(RR_counts)
    for i,j in enumerate(J_regions):
        print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
        
        # Compute pair counts between data jackknife and random survey
        filt = np.where(dJ==j)
        if len(filt[0])>0:
            cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],
                                   weights1=dW[filt],weight_type='pair_product', X2=rX, Y2=rY, 
                                   Z2 = rZ, weights2 = rW, verbose=False,periodic=True)
            DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
        
        # Compute pair coutnts between random jackknife and data survey
        filt = np.where(rJ==j)
        if len(filt[0])>0:
            cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                   weights1=dW,weight_type='pair_product', X2=rX[filt], Y2=rY[filt], 
                                   Z2 = rZ[filt], weights2 = rW[filt], verbose=False,periodic=True)
            DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
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
    print("Finished after %d seconds"%(time.time()-init))
    
from Corrfunc.utils import convert_3d_counts_to_cf

# Now compute correlation function
xi_function = np.zeros_like(RR_counts)
for j in range(N_jack):
    print("CHECK N_GAL NORMALIZATION IN XI_JACK")
    xi_function[j] = convert_3d_counts_to_cf(N_gal,N_gal,N_rand,N_rand,DD_counts[j],
                                             DR_counts[j],DR_counts[j],RR_counts[j])

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

outname='xi_jack_n%d_m%d_j%d_11.dat'%(nrbins,nmu_bins,N_jack)
print("Saving correlation function")
with open(outdir+outname,"w+") as outfile:
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

print("TESTING - SAVE ALL")
np.savez(outdir+'all_xi_jackknife.npz',RR=RR_counts,DR=DR_counts,DD=DD_counts,xi=xi_function)
print("Saved pair counts to %s"%(outdir+'all_xi_jackknife.npz'))
