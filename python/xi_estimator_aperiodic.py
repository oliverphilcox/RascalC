## Script to compute estimate of the correlation functions xi(r,mu) for an entire survey via the Landay-Szalay estimator for one or two sets of random particles.
## This creates three correlation functions - for field 1 x field 1, field 1 x field 2, field 2 x field 2.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=9:
    if len(sys.argv)!=12:
        if len(sys.argv)!=15:
            if len(sys.argv)!=10:
                print("Usage: python xi_estimator_aperiodic.py {GALAXY_FILE} {RANDOM_FILE_DR} {RANDOM_FILE_RR}  {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {OUTPUT_DIR} [{GALAXY_FILE_2} {RANDOM_FILE_2_DR} {RANDOM_FILE_2_RR}] [{RR_counts_11} {RR_counts_12} {RR_counts_22}]")
                sys.exit()

Dname1 = str(sys.argv[1])
Rname1_DR = str(sys.argv[2])
Rname1_RR = str(sys.argv[3])
binfile = str(sys.argv[4])
mu_max = float(sys.argv[5])
nmu_bins = int(sys.argv[6])
nthreads = int(sys.argv[7])
outdir=str(sys.argv[8])

if len(sys.argv)>=12:
    multifield = True
    Dname2 = str(sys.argv[9])
    Rname2_DR = str(sys.argv[10])
    Rname2_RR = str(sys.argv[11])
else:
    multifield = False

if len(sys.argv)==15:
    print("Using pre-defined RR counts")
    RRname11=str(sys.argv[12])
    RRname12=str(sys.argv[13])
    RRname22=str(sys.argv[14])
else:
    if len(sys.argv)==10:
        RRname11 = str(sys.argv[9])
    else:
        RRname11=""
    RRname12=""
    RRname22=""

## First read in weights and positions:
dtype = np.double

def reader(filename,only_length=False):
    """Read in input file"""
    print("Counting lines in file %s"%filename)
    total_lines=0
    for n,line in enumerate(open(filename,"r")):
        total_lines+=1

    if only_length:
        return total_lines

    X,Y,Z,W=[np.zeros(total_lines) for _ in range(4)]

    for n, line in enumerate(open(filename, 'r')):
        if n%1000000==0:
            print("Reading line %d of %d from file %s" %(n,total_lines,filename))
        split_line=np.array(line.split(), dtype=float)
        X[n]=split_line[0];
        Y[n]=split_line[1];
        Z[n]=split_line[2];
        W[n]=split_line[3];

    return X,Y,Z,W

## Read galaxy files
data1 = reader(Dname1)

## Check for periodicity
xrange = max(data1[0])-min(data1[0])
yrange = max(data1[1])-min(data1[1])
zrange = max(data1[2])-min(data1[2])

if (np.abs(xrange-yrange)/xrange<1e-2) and (np.abs(zrange-xrange)/xrange<1e-2):
    print('Data set 1 seems to be periodic! The xi_estimator_periodic function should be used here!')
    sys.exit()

if multifield:
    data2 = reader(Dname2)

    ## Check for periodicity
    xrange = max(data2[0])-min(data2[0])
    yrange = max(data2[1])-min(data2[1])
    zrange = max(data2[2])-min(data2[2])

    if (np.abs(xrange-yrange)/xrange<1e-2) and (np.abs(zrange-xrange)/xrange<1e-2):
        print('Data set 2 seems to be periodic! The xi_estimator_periodic function should be used here!')
        sys.exit()

## Read DR random files
random1_DR = reader(Rname1_DR)

# Total particle numbers
N_rand1_DR = len(random1_DR[0])
N_gal1 = len(data1[0])

if multifield:
    random2_DR = reader(Rname2_DR)

    # Total particle numbers
    N_rand2_DR = len(random2_DR[0])
    N_gal2 = len(data2[0])

# Read RR random files
if multifield:
    # Read RR random files
    if len(RRname11)==0:
        # if RR counts are not provided
        if Rname1_DR!=Rname1_RR:
            random1_RR = reader(Rname1_RR)
        else:
            random1_RR = random1_DR
        if Rname2_DR!=Rname2_RR:
            random2_RR = reader(Rname2_RR)
        else:
            random2_RR = random2_DR
        N_rand1_RR = len(random1_RR[0])
        N_rand2_RR = len(random2_RR[0])
    else:
        # empty placeholders
        random1_RR = []
        random2_RR = []
        N_rand1_RR = reader(Rname1_RR,only_length=True)
        N_rand2_RR = reader(Rname2_RR,only_length=True)
else:
    if len(RRname11)==0:
        # if RR counts are not provided
        if Rname1_DR!=Rname1_RR:
            # if not already read in
            random1_RR = reader(Rname1_RR)
        else:
            random1_RR = random1_DR
        N_rand1_RR = len(random1_RR[0])
    else:
        # empty placeholders
        random1_RR = []
        N_rand1_RR = reader(Rname1_RR,only_length=True)

print("Number of random particles in field 1: %.1e (DR) %.1e (RR)"%(N_rand1_DR,N_rand1_RR))
print("Number of galaxies particles in field 1: %.1e"%N_gal1)

if multifield:
    print("Number of random particles in field 2: %.1e (DR) %.1e (RR)"%(N_rand2_DR,N_rand2_RR))
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

def compute_xi(random1_RR,random1_DR,data1,N_gal,N_rand_RR, N_rand_DR,random2_RR=None,random2_DR=None,data2=None,N_gal2=None,N_rand2_RR=None,N_rand2_DR=None,cross_term=False,RRname="",verbose=False):
    """ Compute the Xi estimate for a given pair of fields. If cross_term=True, we use both fields. If RRname is non-empty we use this instead of recomputing RR pair counts."""

    # Read in fields
    rX_DR,rY_DR,rZ_DR,rW_DR = random1_DR
    if len(random1_RR)>0:
        # if these files have been read in
        rX_RR,rY_RR,rZ_RR,rW_RR = random1_RR
    dX,dY,dZ,dW = data1
    if cross_term:
        rX2_DR,rY2_DR,rZ2_DR,rW2_DR = random2_DR
        if len(random2_RR)>0:
            rX2_RR,rY2_RR,rZ2_RR,rW2_RR = random2_RR
        dX2,dY2,dZ2,dW2 = data2

    import time
    init = time.time()

    # Compute RR, DR and DD counts for the non-periodic case (measuring mu from the radial direction)
    print("Using non-periodic input data");

    # Convert coordinates to spherical coordinates
    r_com_dist_DR,r_Ra_DR,r_Dec_DR = coord_transform(rX_DR,rY_DR,rZ_DR);
    d_com_dist,d_Ra,d_Dec = coord_transform(dX,dY,dZ);
    if cross_term:
        r_com_dist2_DR,r_Ra2_DR,r_Dec2_DR = coord_transform(rX2_DR,rY2_DR,rZ2_DR);
        d_com_dist2,d_Ra2,d_Dec2 = coord_transform(dX2,dY2,dZ2);

    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks

    # Now compute RR counts
    if len(RRname)!=0:
        RR_counts = np.loadtxt(RRname) # read pre-computed RR counts
        RR_counts/=(N_rand_RR-1.)*N_rand_RR
        if len(RR_counts)!=nrbins*nmu_bins:
            raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts.")
    else:
        r_com_dist_RR,r_Ra_RR,r_Dec_RR = coord_transform(rX_RR,rY_RR,rZ_RR);
        print("Computing RR pair counts")
        if cross_term:
            r_com_dist2_RR,r_Ra2_RR,r_Dec2_RR = coord_transform(rX2_RR,rY2_RR,rZ2_RR);
            tmpRR=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,
                             RA2=r_Ra2_RR,DEC2 = r_Dec2_RR, CZ2 = r_com_dist2_RR, weights2 = rW2_RR, weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
            RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg']/ (np.sum(rW_RR)*np.sum(rW2_RR))
        else:
            tmpRR=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
            RR_counts = tmpRR[:]['npairs']*tmpRR[:]['weightavg']/ (np.sum(rW_RR)**2.)
    print("Finished after %d seconds"%(time.time()-init))

    # Now compute DR counts
    print("Computing DR pair counts")
    if cross_term:
        tmpD1R2 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',
                            RA2=r_Ra2_DR, DEC2=r_Dec2_DR, CZ2 = r_com_dist2_DR, weights2 = rW2_DR, verbose=verbose,is_comoving_dist=True)
        tmpD2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2,d_Dec2,d_com_dist2,weights1=dW2,weight_type='pair_product',
                            RA2=r_Ra_DR, DEC2=r_Dec_DR, CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=verbose,is_comoving_dist=True)
        D1R2_counts = tmpD1R2[:]['npairs']*tmpD1R2[:]['weightavg']/(np.sum(rW2_DR)*np.sum(dW))
        D2R1_counts = tmpD2R1[:]['npairs']*tmpD2R1[:]['weightavg']/(np.sum(dW2)*np.sum(rW_DR))
    else:
        tmpDR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',
                            RA2=r_Ra_DR, DEC2=r_Dec_DR, CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=verbose,is_comoving_dist=True)
        DR_counts = tmpDR[:]['npairs']*tmpDR[:]['weightavg']/(np.sum(rW_DR)*np.sum(dW))
    print("Finished after %d seconds"%(time.time()-init))

    # Now compute DD counts
    print("Compute DD pair counts")
    if cross_term:
        tmpDD=DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                          RA2=d_Ra2, DEC2=d_Dec2, CZ2=d_com_dist2, weights2=dW2, weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
        DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']/(np.sum(dW)*np.sum(dW2))
    else:
        tmpDD=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,weight_type='pair_product',verbose=verbose,is_comoving_dist=True)
        DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']/(np.sum(dW)**2.)
    print("Finished after %d seconds"%(time.time()-init))

    # Now use Landay-Szalay estimator:
    if cross_term:
        xi_function = DD_counts/RR_counts - D1R2_counts/RR_counts - D2R1_counts/RR_counts + 1.
    else:
        xi_function = DD_counts/RR_counts - 2.*DR_counts/RR_counts + 1.

    return xi_function.reshape(nrbins,nmu_bins)

## Now compute correlation functions.

print("\nCOMPUTING xi_11 CORRELATION\n")
xi_11=compute_xi(random1_RR,random1_DR,data1,N_gal1,N_rand1_RR,N_rand1_DR,cross_term=False,RRname=RRname11,verbose=False)

if multifield:
    print("\nCOMPUTING xi_12 CORRELATION\n")
    xi_12=compute_xi(random1_RR,random1_DR,data1,N_gal1,N_rand1_RR,N_rand1_DR,random2_RR,random2_DR,data2,N_gal2,N_rand2_RR, N_rand2_DR,cross_term=True,RRname=RRname12,verbose=False)
    print("\nCOMPUTING xi_22 CORRELATION\n")
    xi_22=compute_xi(random2_RR,random2_DR,data2,N_gal2,N_rand2_RR, N_rand2_DR,cross_term=False,RRname=RRname22,verbose=False)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

## Now save to file
if multifield:
    suffices = ['11','12','22']
    xi_files = [xi_11,xi_12,xi_22]
else:
    suffices = ['11']
    xi_files = [xi_11]

import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

for index in range(len(suffices)):
    outname='xi_n%d_m%d_%s.dat'%(nrbins,nmu_bins,suffices[index])
    print("Saving %s correlation function to %s"%(suffices[index],outdir+outname))
    with open(outdir+outname,"w+") as outfile:
        for r in mean_bins:
            outfile.write("%.8e "%r)
        outfile.write("\n")
        for mu in mean_mus:
            outfile.write("%.8e "%mu)
        outfile.write("\n")
        for i in range(nrbins):
            for j in range(nmu_bins):
                outfile.write("%.8e "%xi_files[index][i,j])
            outfile.write("\n")

print("All correlation functions computed successfully.")

print("NB: Number of galaxies in field 1 is %d"%N_gal1)

if multifield:
    print("NB: Number of galaxies in field 2 is %d"%N_gal2)
