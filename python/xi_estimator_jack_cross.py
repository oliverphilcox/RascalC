## Script to compute jackknife estimates of the correlation function xi^J(r,mu) via the Landay-Szelay estimator for a single set of random particles.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.
## This must be binned in the same binning as the desired covariance matrix

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=11:
    if len(sys.argv)!=14:
        print("Usage: python xi_estimator_jack_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1} {RANDOM_FILE_2} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts_11} {RR_counts_12} {RR_counts_22}]")
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
    RRname11=str(sys.argv[11])  
    RRname12=str(sys.argv[12])
    RRname22=str(sys.argv[13])
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
    
    X,Y,Z,W,J=[np.zeros(total_lines) for _ in range(5)]
    
    for n, line in enumerate(open(filename, 'r')):
        if n%1000000==0:
            print("Reading line %d of %d from file %s" %(n,total_lines,filename))
        split_line=np.array(line.split(" "), dtype=float) 
        X[n]=split_line[0];
        Y[n]=split_line[1];
        Z[n]=split_line[2];
        W[n]=split_line[3];
        J[n]=int(split_line[4]);
    return X,Y,Z,W,J

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


# Determine number of jackknifes
J_regions = np.unique(np.concatenate([random1[4],random2[4],data1[4],data2[4]]))
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
print('%s radial bins are used in this file in the range [%d,%d].' %(nrbins,all_bins[0,0],all_bins[-1,1])

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

def compute_xi(random1,data1,random2=None,data2=None,cross_term=False,RRname="",verbose=False):
    """ Compute the Xi jackknife estimates for a given pair of fields. If cross_term=True, we use both fields. If RRname is non-empty we use this instead of recomputing RR pair counts. Estimates are NaN if there's no particles in the jackknife."""
        
    # Read in fields
    rX,rY,rZ,rW,rJ = random1
    dX,dY,dZ,dW,dJ = data1
    if cross_term:
        rX2,rY2,rZ2,rW2,rJ2 = random2
        dX2,dY2,dZ2,dW2,dJ2 = data2
        
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
        RR_counts = np.zeros([N_jack,nrbins*nmu_bins])    
        if len(RRname)!=0:
            # Load pre-existing counts if present
            RRfile = np.loadtxt(RRname) # read pre-computed RR counts
            if len(RRfile[0,:])!=(1+nrbins*nmu_bins):
                raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
            if len(RR_counts[:,0])!=N_jack:
                raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
            for jk in range(N_jack):
                RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
        else:
            # Compute RR pair counts
            for i,j in enumerate(J_regions):
                # Compute pair counts between jackknife region and entire survey regions
                print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
                if cross_term:
                    filt = np.where(rJ==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra2,r_Dec2,r_com_dist2,weights1=rW2,
                                                weight_type='pair_product',RA2=r_Ra[filt],DEC2=r_Dec[filt],CZ2=r_com_dist[filt],
                                                weights2=rW2[filt],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    filt = np.where(rJ2==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra,r_Dec,r_com_dist,weights1=rW,
                                           weight_type='pair_product',RA2=r_Ra2[filt],DEC2=r_Dec2[filt],CZ2=r_com_dist2[filt],
                                           weights2=rW2[filt],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                else:
                    filt = np.where(rJ==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra,r_Dec,r_com_dist,weights1=rW,
                                                weight_type='pair_product',RA2=r_Ra[filt],DEC2=r_Dec[filt],CZ2=r_com_dist[filt],
                                                weights2=rW[filt],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                
        print("Finished RR pair counts after %d seconds"%(time.time()-init))
        
        if cross_term:
            D1R2_counts = np.zeros_like(RR_counts)
            D2R1_counts = np.zeros_like(RR_counts)
        else:
            DR_counts = np.zeros_like(RR_counts)
        
        # Now compute DR counts
        for i,j in enumerate(J_regions):
            print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if cross_term:
                
                ## Compute D1R2 term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_D1R2 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],
                                             weights1=dW[filt],weight_type='pair_product',RA2=r_Ra2,DEC2=r_Dec2,
                                             CZ2 = r_com_dist2, weights2 = rW2, verbose=verbose,is_comoving_dist=True)
                    D1R2_counts[i,:] +=cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                filt = np.where(rJ2==j)
                if len(filt[0])>0:
                    cross_D1R2 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                             weights1=dW,weight_type='pair_product',RA2=r_Ra2[filt],DEC2=r_Dec2[filt],
                                             CZ2 = r_com_dist2[filt], weights2 = rW2[filt], verbose=verbose,is_comoving_dist=True)
                    D1R2_counts[i,:] +=cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                
                ## Compute D2R1 term
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2[filt],d_Dec2[filt],d_com_dist2[filt],
                                             weights1=dW2[filt],weight_type='pair_product',RA2=r_Ra,DEC2=r_Dec,
                                             CZ2 = r_com_dist, weights2 = rW, verbose=verbose,is_comoving_dist=True)
                    D2R1_counts[i,:] +=cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                filt = np.where(rJ==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2,d_Dec2,d_com_dist2,
                                             weights1=dW2,weight_type='pair_product',RA2=r_Ra[filt],DEC2=r_Dec[filt],
                                             CZ2 = r_com_dist[filt], weights2 = rW[filt], verbose=verbose,is_comoving_dist=True)
                    D2R1_counts[i,:] +=cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                
            else:
                
                # Just compute DR term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],
                                        weights1=dW[filt],weight_type='pair_product', RA2=r_Ra, DEC2=r_Dec, 
                                        CZ2 = r_com_dist, weights2 = rW, verbose=verbose,is_comoving_dist=True)
                    DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                filt = np.where(rJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                        weights1=dW,weight_type='pair_product', RA2=r_Ra[filt], DEC2=r_Dec[filt], 
                                        CZ2 = r_com_dist[filt], weights2 = rW[filt], verbose=verbose,is_comoving_dist=True)
                    DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                    
        print("Finished DR pair counts after %d seconds"%(time.time()-init))
            
        # Now compute DD counts
        DD_counts = np.zeros_like(RR_counts)
        for i,j in enumerate(J_regions):
            if cross_term:
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                                            weight_type='pair_product',RA2=d_Ra2[filt],DEC2=d_Dec2[filt],CZ2=d_com_dist2[filt],
                                            weights2=dW2[filt],verbose=verbose,is_comoving_dist=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],weights1=dW[filt],
                                            weight_type='pair_product',RA2=d_Ra2,DEC2=d_Dec2,CZ2=d_com_dist2,
                                            weights2=dW2,verbose=verbose,is_comoving_dist=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
            else:
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                                            weight_type='pair_product',RA2=d_Ra[filt],DEC2=d_Dec[filt],CZ2=d_com_dist[filt],
                                            weights2=dW[filt],verbose=verbose,is_comoving_dist=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                
        print("Finished DD pair counts after %d seconds"%(time.time()-init))
    
    else:
        
        # Compute xi for the periodic case (measuring mu from the Z-axis)
        print("Using periodic input data")
        from Corrfunc.mocks.DDsmu import DDsmu
        
        import time
        init = time.time()
        
        # Now compute RR counts
        RR_counts = np.zeros([N_jack,nrbins*nmu_bins])    
        if len(RRname)!=0:
            # Load pre-existing counts if present
            RRfile = np.loadtxt(RRname) # read pre-computed RR counts
            if len(RRfile[0,:])!=(1+nrbins*nmu_bins):
                raise Exception("Incorrect number of bins in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
            if len(RR_counts[:,0])!=N_jack:
                raise Exception("Incorrect number of jackknives in RR file. Either provide the relevant file or recompute RR pair counts for each unrestricted jackknife.")
            for jk in range(N_jack):
                RR_counts[jk,:] = RRfile[jk,1:] # first index is jackknife number usually
        else:
            # Compute RR pair counts
            for i,j in enumerate(J_regions):
                # Compute pair counts between jackknife region and entire survey regions
                print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
                if cross_term:
                    filt = np.where(rJ==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX2,rY2,rZ2,weights1=rW2,
                                                weight_type='pair_product',X2=rX[filt],Y2=rY[filt],Z2=rZ[filt],
                                                weights2=rW2[filt],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    filt = np.where(rJ2==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX,rY,rZ,weights1=rW,
                                           weight_type='pair_product',X2=rX2[filt],Y2=rY2[filt],Z2=rZ2[filt],
                                           weights2=rW2[filt],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                else:
                    filt = np.where(rJ==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX,rY,rZ,weights1=rW,
                                                weight_type='pair_product',X2=rX[filt],Y2=rY[filt],Z2=rZ[filt],
                                                weights2=rW[filt],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
        print("Finished RR pair counts after %d seconds"%(time.time()-init))
        
        if cross_term:
            D1R2_counts = np.zeros_like(RR_counts)
            D2R1_counts = np.zeros_like(RR_counts)
        else:
            DR_counts = np.zeros_like(RR_counts)
        
        # Now compute DR counts
        for i,j in enumerate(J_regions):
            print("Computing DR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
            if cross_term:
                
                ## Compute D1R2 term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_D1R2 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],
                                             weights1=dW[filt],weight_type='pair_product',X2=rX2,Y2=rY2,
                                             Z2 = rZ2, weights2 = rW2, verbose=verbose,periodic=True)
                    D1R2_counts[i,:] +=cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                filt = np.where(rJ2==j)
                if len(filt[0])>0:
                    cross_D1R2 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                             weights1=dW,weight_type='pair_product',X2=rX2[filt],Y2=rY2[filt],
                                             Z2 = rZ2[filt], weights2 = rW2[filt], verbose=verbose,periodic=True)
                    D1R2_counts[i,:] +=cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                
                ## Compute D2R1 term
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX2[filt],dY2[filt],dZ2[filt],
                                             weights1=dW2[filt],weight_type='pair_product',X2=rX,Y2=rY,
                                             Z2 = rZ, weights2 = rW, verbose=verbose,periodic=True)
                    D2R1_counts[i,:] +=cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                filt = np.where(rJ==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX2,dY2,dZ2,
                                             weights1=dW2,weight_type='pair_product',X2=rX[filt],Y2=rY[filt],
                                             Z2 = rZ[filt], weights2 = rW[filt], verbose=verbose,periodic=True)
                    D2R1_counts[i,:] +=cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                
            else:
                
                # Just compute DR term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],
                                        weights1=dW[filt],weight_type='pair_product', X2=rX, Y2=rY, 
                                        Z2 = rZ, weights2 = rW, verbose=verbose,periodic=True)
                    DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                filt = np.where(rJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                        weights1=dW,weight_type='pair_product', X2=rX[filt], Y2=rY[filt], 
                                        Z2 = rZ[filt], weights2 = rW[filt], verbose=verbose,periodic=True)
                    DR_counts[i,:] += cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                    
        print("Finished DR pair counts after %d seconds"%(time.time()-init))
            
        # Now compute DD counts
        DD_counts = np.zeros_like(RR_counts)
        for i,j in enumerate(J_regions):
            if cross_term:
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,
                                            weight_type='pair_product',X2=dX2[filt],Y2=dY2[filt],Z2=dZ2[filt],
                                            weights2=dW2[filt],verbose=verbose,periodic=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],weights1=dW[filt],
                                            weight_type='pair_product',X2=dX2,Y2=dY2,Z2=dZ2,
                                            weights2=dW2,verbose=verbose,periodic=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
            else:
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,
                                            weight_type='pair_product',X2=dX[filt],Y2=dY[filt],Z2=dZ[filt],
                                            weights2=dW[filt],verbose=verbose,periodic=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                
        print("Finished DD pair counts after %d seconds"%(time.time()-init))
    
    # Compute normalizations
    N_RR = np.sum(RR_counts,axis=1)
    N_DD = np.sum(DD_counts,axis=1)
        
    ## Now compute correlation function
    if cross_term:
        N_D1R2 = np.sum(D1R2_counts,axis=1)
        N_D2R1 = np.sum(D2R1_counts,axis=1)
        
        # Now use Landay-Szelay estimator:
        for j in range(N_jack):
            xi_function[j] = DD_counts[j]/RR_counts[j]*N_RR[j]/N_DD[j] - D1R2_counts[j]/RR_counts[j]*N_RR[j]/N_D1R2[j] - D2R1_counts[j]/RR_counts[j]*N_RR[j]/N_D2R1[j] + 1.
    else:
        N_DR = np.sum(DR_counts,axis=1)
        
        # Now use Landay-Szelay estimator:
        for j in range(N_jack):
            xi_function[j] = DD_counts[j]/RR_counts[j]*N_RR[j]/N_DD[j] - 2.*DR_counts[j]/RR_counts[j]*N_RR[j]/N_DR[j] + 1.

    return xi_function

## Now compute correlation functions.

print("\nCOMPUTING xi_11 CORRELATION\n")
xi_11=compute_xi(random1,data1,cross_term=False,RRname=RRname11,verbose=False)
print("\nCOMPUTING xi_12 CORRELATION\n")
xi_12=compute_xi(random1,data1,random2,data2,cross_term=True,RRname=RRname12,verbose=False)
print("\nCOMPUTING xi_22 CORRELATION\n")
xi_22=compute_xi(random2,data2,cross_term=False,RRname=RRname22,verbose=False)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

## Now save to file
suffices = ['11','12','22']
xi_files = [xi_11,xi_12,xi_22]

import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

for index in range(3):
    outname='xi_jack_n%d_m%d_%s.dat'%(nrbins,nmu_bins,suffices[index])
    print("Saving %s correlation function to %s"%(suffices[index],outdir+outname))
    with open(outdir+outname,"w+") as outfile:
        for r in mean_bins:
            outfile.write("%.8e "%r)
        outfile.write("\n")
        for mu in mean_mus:
            outfile.write("%.8e "%mu)
        outfile.write("\n")
        for i in range(N_jack):
            for j in range(nrbins*nmu_bins):
                outfile.write("%.8e "%xi_files[index][i,j])
            outfile.write("\n")

print("All correlation functions computed successfully.")
