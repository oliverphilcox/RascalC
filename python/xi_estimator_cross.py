## Script to compute estimate of the correlation functions xi(r,mu) for an entire survey via the Landay-Szelay estimator for two sets of random particles.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.
## This creates three correlation functions - for field 1 x field 1, field 1 x field 2, field 2 x field 2.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=11:
    if len(sys.argv)!=14:
        print("Usage: python xi_estimator_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1} {RANDOM_FILE_2} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts_11} {RR_counts_12} {RR_counts_22}]")
        sys.exit()
Dname1 = str(sys.argv[1])
Dname2 = str(sys.argv[2])
Rname1 = str(sys.argv[3])
Rname2 = str(sys.argv[4])
binfile = str(sys.argv[5])
mu_max = float(sys.argv[6])
nmu_bins = int(sys.argv[7])
nthreads = int(sys.argv[8])
periodic = int(sys.argv[9])
outdir=str(sys.argv[10])

if len(sys.argv)==14:
    print("Using pre-defined RR counts")
    RRname11=str(sys.argv[9])
    RRname12=str(sys.argv[10])
    RRname22=str(sys.argv[11])
else:
    RRname11=""
    RRname12=""
    RRname22=""

## First read in weights and positions:
dtype = np.double 

def reader(filename):
    """Read in input file"""
    print("Counting lines in file %s"%filename)
    total_lines=0
    for n,line in enumerate(open(filename,"r")):
        total_lines+=1
    
    X,Y,Z,W=[np.zeros(total_lines) for _ in range(4)]
    
    for n, line in enumerate(open(filename, 'r')):
        if n%1000000==0:
            print("Reading line %d of %d from file %s" %(n,total_lines,filename))
        split_line=np.array(line.split(" "), dtype=float) 
        X[n]=split_line[0];
        Y[n]=split_line[1];
        Z[n]=split_line[2];
        W[n]=split_line[3];

    return X,Y,Z,W

## Read random files
random1 = reader(Rname1)
random2 = reader(Rname2)

## Read galaxy files
data1 = reader(Dname1)
data2 = reader(Dname2)

# Total particle numbers
N_rand1 = len(random1[0])
N_rand2 = len(random2[0]) 
N_gal1 = len(data1[0])
N_gal2 = len(data2[0])

print("Number of random particles in field 1: %.1e"%N_rand1)
print("Number of galaxies particles in field 1: %.1e"%N_gal1)
print("Number of random particles in field 2: %.1e"%N_rand2)
print("Number of galaxies particles in field 2: %.1e"%N_gal2)

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

## Coordinate transformations
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

def compute_xi(random1,data1,N_gal,N_rand,random2=None,data2=None,N_gal2=None,N_rand2=None,cross_term=False,RRname="",verbose=False):
    """ Compute the Xi estimate for a given pair of fields. If cross_term=True, we use both fields. If RRname is non-empty we use this instead of recomputing RR pair counts."""
        
    # Read in fields
    rX,rY,rZ,rW = random1
    dX,dY,dZ,dW = data1
    if cross_term:
        rX2,rY2,rZ2,rW2 = random2
        dX2,dY2,dZ2,dW2 = data2
        
    import time
    init = time.time()
        
    if not periodic:
        # Compute RR, DR and DD counts for the non-periodic case (measuring mu from the radial direction)
        print("Using non-periodic input data");
        
        # Convert coordinates to spherical coordinates
        r_com_dist,r_Ra,r_Dec = coord_transform(rX,rY,rZ);
        d_com_dist,d_Ra,d_Dec = coord_transform(dX,dY,dZ);
        if cross_term:
            r_com_dist2,r_Ra2,r_Dec2 = coord_transform(rX2,rY2,rZ2);
            d_com_dist2,d_Ra2,d_Dec2 = coord_transform(dX2,dY2,dZ2);

        from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
        
        # Now compute RR counts
        if len(RRname)!=0:
            RR_counts = np.loadtxt(RRname) # read pre-computed RR counts
            if len(RR_counts)!=nrbins*nmu_bins:
                raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts.")
        else:
            print("Computing RR pair counts")
            if cross_term:
               tmpRR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra,r_Dec,com_dist,weights1=rW,
                                 RA2=r_Ra2,DEC2 = r_Dec2, CZ2 = r_com_dist2, weights2 = rW2, weight_type='pair_product',verbose=verbose,is_comoving_dist=True) 
            else:
                tmpRR=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,r_Ra,r_Dec,r_com_dist,weights1=rW,weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
            RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg']
        print("Finished after %d seconds"%(time.time()-init))
        
        # Now compute DR counts
        print("Computing DR pair counts")
        if cross_term:
            tmpD1R2 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',
                                RA2=r_Ra2, DEC2=r_Dec2, CZ2 = r_com_dist2, weights2 = rW2, verbose=verbose,is_comoving_dist=True)
            tmpD2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2,d_Dec2,d_com_dist2,weights1=dW2,weight_type='pair_product',
                                RA2=r_Ra, DEC2=r_Dec, CZ2 = r_com_dist, weights2 = rW, verbose=verbose,is_comoving_dist=True)
            D1R2_counts = tmpD1R2[:]['npairs']*tmpD1R2[:]['weightavg']
            D2R1_counts = tmpD2R1[:]['npairs']*tmpD2R1[:]['weightavg']
        else:
            tmpDR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',
                                RA2=r_Ra, DEC2=r_Dec, CZ2 = r_com_dist, weights2 = rW, verbose=verbose,is_comoving_dist=True)
            DR_counts = tmpDR[:]['npairs']*tmpDR[:]['weightavg']
        print("Finished after %d seconds"%(time.time()-init))
        
        # Now compute DD counts
        print("Compute DD pair counts")
        if cross_term:
            tmpDD=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                              RA2=d_Ra2, DEC2=d_Dec2, CZ2=d_com_dist2, weights2=dW2, weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
        else:
            tmpDD=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
        DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']
        print("Finished after %d seconds"%(time.time()-init))
        
    else:
        # Compute RR counts for the periodic case (measuring mu from the Z-axis)
        print("Using periodic input data");
        from Corrfunc.theory.DDsmu import DDsmu
        
        # Now compute RR counts
        if len(RRname)!=0:
            RR_counts = np.loadtxt(RRname) # read pre-computed RR counts
            if len(RR_counts)!=nrbins*nmu_bins:
                raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts.")
        else:
            print("Computing RR pair counts")
            if cross_term:
                tmpRR=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX,rY,rZ,weights1=rW,
                            X2=rX2,Y2=rY2,Z2=rZ2,W2=rW2, verbose=verbose,periodic=True)
            else:
                tmpRR=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,rX,rY,rZ,weights1=rW,weight_type='pair_product',verbose=verbose,periodic=True)
            RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg']
            print("Finished after %d seconds"%(time.time()-init))
            
        # Now compute DR counts
        print("Computing DR pair counts")
        if cross_term:
            tmpD1R2 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',
                            X2=rX2, Y2=rY2, Z2 = rZ2, weights2 = rW2, verbose=verbose,periodic=True)
            tmpD2R1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX2,dY2,dZ2,weights1=dW2,weight_type='pair_product',
                            X2=rX, Y2=rY, Z2 = rZ, weights2 = rW, verbose=verbose,periodic=True)
            D1R2_counts = tmpD1R2[:]['npairs']*tmpD1R2[:]['weightavg']
            D2R1_counts = tmpD2R1[:]['npairs']*tmpD2R1[:]['weightavg']
        else:
            tmpDR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',
                            X2=rX, Y2=rY, Z2 = rZ, weights2 = rW, verbose=verbose,periodic=True)
            DR_counts = tmpDR[:]['npairs']*tmpDR[:]['weightavg']
        print("Finished after %d seconds"%(time.time()-init))
        
        # Now compute DD counts
        print("Compute DD pair counts")
        if cross_term:
            tmpDD=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,
                        X2=dX2,Y2=dY2,Z2=dZ2,weights2=dW2,weight_type='pair_product',verbose=verbose,periodic=True)
        DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']
        print("Finished after %d seconds"%(time.time()-init))
        
    # Now reconstruct xi function:
    from Corrfunc.utils import convert_3d_counts_to_cf

    # Now compute correlation function
    if cross_term:
        xi_function = convert_3d_counts_to_cf(N_gal,N_gal2,N_rand,N_rand2,DD_counts,D1R2_counts,D2R1_counts,RR_counts)
    else:
        xi_function = convert_3d_counts_to_cf(N_gal,N_gal,n_rand,N_rand,DD_counts,DR_counts,DR_counts,RR_counts)
        
    return xi_function.reshape(nrbins,nmu_bins)


## Now compute correlation functions.

print("\nCOMPUTING xi_11 CORRELATION\n")
xi_11=compute_xi(random1,data1,N_gal1,N_rand1,cross_term=False,RRname=RRname11,verbose=False)
print("\nCOMPUTING xi_12 CORRELATION\n")
xi_12=compute_xi(random1,data1,N_gal1,N_rand1,random2,data2,N_gal2,N_rand2,cross_term=True,RRname=RRname12,verbose=False)
print("\nCOMPUTING xi_22 CORRELATION\n")
xi_22=compute_xi(random2,data2,N_gal2,N_rand2,cross_term=False,RRname=RRname22,verbose=False)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

## Now save to file
suffices = ['11','12','22']
xi_files = [xi_11,xi_12,xi_22]

import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

for index in range(3):
    outname='xi_n%d_m%d_%s.dat'%(nrbins,nmu_bins,suffices[index])
    print("Saving %d-th correlation function to %s"%(i+1,outdir+outname))
    with open(outdir+outname,"w+") as outfile:
        for r in mean_bins:
            outfile.write("%.8e "%r)
        outfile.write("\n")
        for mu in mean_mus:
            outfile.write("%.8e "%mu)
        for i in range(nrbins):
            for j in range(nmu_bins):
                outfile.write("%.8e "%xi_files[index][i,j])
            outfile.write("\n")
            
print("All correlation functions computed successfully.")
            
