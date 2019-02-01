## Script to compute estimate of the single correlation function xi(r,mu) for an entire survey via the Landay-Szelay estimator for a single set of random particles. 
## We allow for two sets of random particles to be used, to allow for different number of random particles in the RR and DR estimation.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=10:
    if len(sys.argv)!=11:
        print("Usage: python xi_estimator.py {GALAXY_FILE} {RANDOM_FILE_DR} {RANDOM_FILE_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts}]")
        sys.exit()
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
dtype = np.double 

# Read first set of randoms
print("Counting lines in DR random file")
total_lines=0
for n, line in enumerate(open(RnameDR, 'r')):
    total_lines+=1

rX_DR,rY_DR,rZ_DR,rW_DR=[np.zeros(total_lines) for _ in range(4)]

print("Reading in DR random data");
for n, line in enumerate(open(RnameDR, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    rX_DR[n]=split_line[0];
    rY_DR[n]=split_line[1];
    rZ_DR[n]=split_line[2];
    rW_DR[n]=split_line[3];
    
N_randDR = len(rX_DR) # number of particles

if len(RRname)==0:
    # Read in RR random file    
    if RnameDR!=RnameRR:
        print("Counting lines in RR random file")
        total_lines=0
        for n, line in enumerate(open(RnameRR, 'r')):
            total_lines+=1

        rX_RR,rY_RR,rZ_RR,rW_RR=[np.zeros(total_lines) for _ in range(4)]

        print("Reading in RR random data");
        for n, line in enumerate(open(RnameRR, 'r')):
            if n%1000000==0:
                print("Reading line %d of %d" %(n,total_lines))
            split_line=np.array(line.split(" "), dtype=float) 
            rX_RR[n]=split_line[0];
            rY_RR[n]=split_line[1];
            rZ_RR[n]=split_line[2];
            rW_RR[n]=split_line[3];

        N_randRR = len(rX_RR) 
    else:
        rX_RR=rX_DR
        rY_RR=rY_DR
        rZ_RR=rZ_DR
        rW_RR=rW_DR
        N_randRR=N_randDR
else:
    print("Counting lines in RR random file")
    N_randRR=0
    for n, line in enumerate(open(RnameRR, 'r')):
        N_randRR+=1
    # empty placeholders
    rX_RR,rY_RR,rZ_RR,rW_RR=[],[],[],[]

print("Counting lines in galaxy file")
total_lines=0
for n, line in enumerate(open(Dname, 'r')):
    total_lines+=1

dX,dY,dZ,dW=[np.zeros(total_lines) for _ in range(4)]

print("Reading in galaxy data");
for n, line in enumerate(open(Dname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float) 
    dX[n]=split_line[0];
    dY[n]=split_line[1];
    dZ[n]=split_line[2];
    dW[n]=split_line[3];

N_gal = len(dX) # number of particles

print("Number of random particles %.1e (DR) %.1e (RR)"%(N_randDR,N_randRR))
print("Number of galaxy particles %.1e"%N_gal)

## Determine number of radial bins in binning file:
print("Counting lines in binfile");
with open(binfile) as f:
    for i, l in enumerate(f):
        pass
nrbins = i + 1
all_bins = np.loadtxt(binfile)
mean_bins=0.5*(all_bins[:,0]+all_bins[:,1])
if all_bins[0,0]>2:
    raise Exception("Radial binfile should extend close to zero")
print('%s radial bins are used in this file.' %nrbins)

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
    if len(RRname)!=0:
        RR_counts = np.loadtxt(RRname) # read pre-computed RR counts
        # these are already normalized
        if len(RR_counts)!=nrbins*nmu_bins:
            raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts.")
    else:
        r_com_dist_RR,r_Ra_RR,r_Dec_RR = coord_transform(rX_RR,rY_RR,rZ_RR);
        print("Computing RR pair counts")
        tmpRR=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,weight_type='pair_product',verbose=False,is_comoving_dist=True)
        RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg'] # sum of weights over bin
        RR_counts/=np.sum(rW_RR)**2.
        print("Finished after %d seconds"%(time.time()-init))
    # Now compute DR counts
    print("Computing DR pair counts")
    tmpDR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',
                        RA2=r_Ra_DR, DEC2=r_Dec_DR, CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=False,is_comoving_dist=True)
    DR_counts = tmpDR[:]['npairs']*tmpDR[:]['weightavg']
    DR_counts/= np.sum(rW_DR)*np.sum(dW)
    print("Finished after %d seconds"%(time.time()-init))
    
    # Now compute DD counts
    print("Compute DD pair counts")
    tmpDD=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',verbose=False,is_comoving_dist=True)
    DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']
    DD_counts/= np.sum(dW)**2.
    print("Finished after %d seconds"%(time.time()-init))
    
else:
    # Compute RR counts for the periodic case (measuring mu from the Z-axis)
    print("Using periodic input data");
    from Corrfunc.theory.DDsmu import DDsmu
    
    import time
    init = time.time()
    
    # Now compute RR counts
    if len(RRname)!=0:
        RR_counts = np.loadtxt(RRname) # read pre-computed RR counts
        # already normalized here
        if len(RR_counts)!=nrbins*nmu_bins:
            raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts.")
    else:
        print("Computing RR pair counts")
        tmpRR=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,rX_RR,rY_RR,rZ_RR,weights1=rW_RR,weight_type='pair_product',verbose=True,periodic=True)
        RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg'] # sum of weights over bin
        RR_counts/= np.sum(rW_RR)**2.
        print("Finished after %d seconds"%(time.time()-init))
    
    # Now compute DR counts
    print("Computing DR pair counts")
    tmpDR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',
                        X2=rX_DR, Y2=rY_DR, Z2 = rZ_DR, weights2 = rW_DR, verbose=True,periodic=True)
    DR_counts = tmpDR[:]['npairs']*tmpDR[:]['weightavg']
    DR_counts/=np.sum(rW_DR)*np.sum(dW)
    print("Finished after %d seconds"%(time.time()-init))
    
    # Now compute DD counts
    print("Compute DD pair counts")
    tmpDD=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',verbose=True,periodic=True)
    DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']
    DD_counts/=np.sum(dW)**2.
    print("Finished after %d seconds"%(time.time()-init))
    
# Now use Landay-Szelay estimator:
xi_function = DD_counts/RR_counts - 2.*DR_counts/RR_counts + 1.

# Now reshape correlation function and save to file.
xi_reshape = xi_function.reshape(nrbins,nmu_bins)

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

outname='xi_n%d_m%d_11.dat'%(nrbins,nmu_bins)
print("Saving correlation function")
with open(outdir+outname,"w+") as outfile:
    for r in mean_bins:
        outfile.write("%.8e "%r)
    outfile.write("\n")
    for mu in mean_mus:
        outfile.write("%.8e "%mu)
    outfile.write("\n");
    for i in range(nrbins):
        for j in range(nmu_bins):
            outfile.write("%.8e "%xi_reshape[i,j])
        outfile.write("\n")
        
print("Correlation function written successfully to %s"%(outdir+outname))

print("NB: Sum of galaxy weights is %.8e"%np.sum(dW))
