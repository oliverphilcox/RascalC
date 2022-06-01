## Script to compute jackknife estimates of the correlation function xi^J(r,mu) via the Landay-Szelay estimator for a single set of random particles.
## If the periodic flag is set, we assume a periodic simulation and measure mu from the Z-axis.
## This must be binned in the same binning as the desired covariance matrix

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=13:
    if len(sys.argv)!=16:
        print("Usage: python xi_estimator_jack_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1_DR} {RANDOM_FILE_1_RR} {RANDOM_FILE_2_DR} {RANDOM_FILE_2_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_jackknife_counts_11} {RR_jackknife_counts_12} {RR_jackknife_counts_22}]")
        sys.exit()
Dname1 = str(sys.argv[1])
Dname2 = str(sys.argv[2])
Rname1_DR = str(sys.argv[3])
Rname1_RR = str(sys.argv[4])
Rname2_DR = str(sys.argv[5])
Rname2_RR = str(sys.argv[6])
binfile = str(sys.argv[7])
mu_max = float(sys.argv[8])
nmu_bins = int(sys.argv[9])
nthreads = int(sys.argv[10])
periodic = int(sys.argv[11])
outdir=str(sys.argv[12])

if len(sys.argv)==16:
    print("Using pre-defined RR counts")
    RRname11=str(sys.argv[13])  
    RRname12=str(sys.argv[14])
    RRname22=str(sys.argv[15])
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

## Read galaxy files
data1 = reader(Dname1)
data2 = reader(Dname2)

## Read DR random files
random1_DR = reader(Rname1_DR)
random2_DR = reader(Rname2_DR)

# Total particle numbers
N_rand1_DR = len(random1_DR[0])
N_rand2_DR = len(random2_DR[0]) 
N_gal1 = len(data1[0])
N_gal2 = len(data2[0])

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

print("Number of random particles in field 1: %.1e (DR) %.1e (RR)"%(N_rand1_DR,N_rand1_RR))
print("Number of galaxies particles in field 1: %.1e"%N_gal1)
print("Number of random particles in field 2: %.1e (DR) %.1e (RR)"%(N_rand2_DR,N_rand2_RR))
print("Number of galaxies particles in field 2: %.1e"%N_gal2)


# Determine number of jackknifes
if len(RRname11)==0:
    J_regions = np.unique(np.concatenate([random1_RR[4],random1_DR[4],random2_RR[4],random2_DR[4],data1[4],data2[4]]))
else:
    J_regions = np.unique(np.concatenate([random1_DR[4],random2_DR[4],data1[4],data2[4]]))
    
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
print('%s radial bins are used in this file in the range [%d,%d].' %(nrbins,all_bins[0,0],all_bins[-1,1]))

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

def compute_xi(random1_RR,random1_DR,data1,N_rand1_RR, N_rand1_DR, N_gal1,random2_RR=None,random2_DR=None,data2=None,N_rand2_RR=None,N_rand2_DR=None,N_gal2=None,cross_term=False,RRname="",verbose=False):
    """ Compute the Xi jackknife estimates for a given pair of fields. If cross_term=True, we use both fields. If RRname is non-empty we use this instead of recomputing RR pair counts. Estimates are NaN if there's no particles in the jackknife."""
        
    # Read in fields
    rX_DR,rY_DR,rZ_DR,rW_DR,rJ_DR = random1_DR
    if len(random1_RR)>0:
        rX_RR,rY_RR,rZ_RR,rW_RR,rJ_RR = random1_RR
    dX,dY,dZ,dW,dJ = data1
    if cross_term:
        rX2_DR,rY2_DR,rZ2_DR,rW2_DR,rJ2_DR = random2_DR
        if len(random2_RR)>0:
            rX2_RR,rY2_RR,rZ2_RR,rW2_RR,rJ2_RR = random2_RR
        dX2,dY2,dZ2,dW2,dJ2 = data2
        
    import time
    init = time.time()
    
    if not periodic:
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
                # NB: these are already normalized 
        else:
            # Compute RR pair counts
            r_com_dist_RR,r_Ra_RR,r_Dec_RR = coord_transform(rX_RR,rY_RR,rZ_RR);
            if cross_term:
                r_com_dist2_RR,r_Ra2_RR,r_Dec2_RR = coord_transform(rX2_RR,rY2_RR,rZ2_RR);
            for i,j in enumerate(J_regions):
                # Compute pair counts between jackknife region and entire survey regions
                print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
                if cross_term:
                    filt = np.where(rJ_RR==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra2_RR,r_Dec2_RR,r_com_dist2_RR,weights1=rW2_RR,
                                                weight_type='pair_product',RA2=r_Ra_RR[filt],DEC2=r_Dec_RR[filt],CZ2=r_com_dist_RR[filt],
                                                weights2=rW_RR[filt],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=0.5*cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    filt2 = np.where(rJ2_RR==j)
                    if len(filt2[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,
                                           weight_type='pair_product',RA2=r_Ra2_RR[filt2],DEC2=r_Dec2_RR[filt2],CZ2=r_com_dist2_RR[filt2],
                                           weights2=rW2_RR[filt2],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=0.5*cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    RR_counts[i,:] /= np.sum(rW2_RR)*np.sum(rW_RR) # normalize by product of sum of weights
                else:
                    filt = np.where(rJ_RR==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,r_Ra_RR,r_Dec_RR,r_com_dist_RR,weights1=rW_RR,
                                                weight_type='pair_product',RA2=r_Ra_RR[filt],DEC2=r_Dec_RR[filt],CZ2=r_com_dist_RR[filt],
                                                weights2=rW_RR[filt],verbose=verbose,is_comoving_dist=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    RR_counts[i,:] /= np.sum(rW_RR)**2 # normalize by product of sum of weights
                
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
                                             weights1=dW[filt],weight_type='pair_product',RA2=r_Ra2_DR,DEC2=r_Dec2_DR,
                                             CZ2 = r_com_dist2_DR, weights2 = rW2_DR, verbose=verbose,is_comoving_dist=True)
                    D1R2_counts[i,:] +=0.5*cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                filt2 = np.where(rJ2_DR==j)
                if len(filt2[0])>0:
                    cross_D1R2 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                             weights1=dW,weight_type='pair_product',RA2=r_Ra2_DR[filt2],DEC2=r_Dec2_DR[filt2],
                                             CZ2 = r_com_dist2_DR[filt2], weights2 = rW2_DR[filt2], verbose=verbose,is_comoving_dist=True)
                    D1R2_counts[i,:] +=0.5*cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                D1R2_counts[i,:] /= np.sum(dW)*np.sum(rW2_DR) # normalize by product of sum of weights
                
                ## Compute D2R1 term
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2[filt],d_Dec2[filt],d_com_dist2[filt],
                                             weights1=dW2[filt],weight_type='pair_product',RA2=r_Ra_DR,DEC2=r_Dec_DR,
                                             CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=verbose,is_comoving_dist=True)
                    D2R1_counts[i,:] +=0.5*cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                filt2 = np.where(rJ_DR==j)
                if len(filt2[0])>0:
                    cross_D2R1 = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra2,d_Dec2,d_com_dist2,
                                             weights1=dW2,weight_type='pair_product',RA2=r_Ra_DR[filt2],DEC2=r_Dec_DR[filt2],
                                             CZ2 = r_com_dist_DR[filt2], weights2 = rW_DR[filt2], verbose=verbose,is_comoving_dist=True)
                    D2R1_counts[i,:] +=0.5*cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                D2R1_counts[i,:] /= np.sum(dW2)*np.sum(rW_DR) # normalize by product of sum of weights
            else:
                # Just compute DR term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt],d_Dec[filt],d_com_dist[filt],
                                        weights1=dW[filt],weight_type='pair_product', RA2=r_Ra_DR, DEC2=r_Dec_DR, 
                                        CZ2 = r_com_dist_DR, weights2 = rW_DR, verbose=verbose,is_comoving_dist=True)
                    DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                filt2 = np.where(rJ_DR==j)
                if len(filt2[0])>0:
                    cross_DR = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,
                                        weights1=dW,weight_type='pair_product', RA2=r_Ra_DR[filt2], DEC2=r_Dec_DR[filt2], 
                                        CZ2 = r_com_dist_DR[filt2], weights2 = rW_DR[filt2], verbose=verbose,is_comoving_dist=True)
                    DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                DR_counts[i,:] /= np.sum(dW)*np.sum(rW_DR) # normalize by product of sum of weights
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
                    DD_counts[i,:]+=0.5*cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                filt2 = np.where(dJ==j)
                if len(filt2[0])>0:
                    cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra[filt2],d_Dec[filt2],d_com_dist[filt2],weights1=dW[filt2],
                                            weight_type='pair_product',RA2=d_Ra2,DEC2=d_Dec2,CZ2=d_com_dist2,
                                            weights2=dW2,verbose=verbose,is_comoving_dist=True)
                    DD_counts[i,:]+=0.5*cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                DD_counts[i,:] /= np.sum(dW)*np.sum(dW2) # normalize by product of sum of weights
            else:
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu_mocks(0,2,nthreads,mu_max,nmu_bins,binfile,d_Ra,d_Dec,d_com_dist,weights1=dW,
                                            weight_type='pair_product',RA2=d_Ra[filt],DEC2=d_Dec[filt],CZ2=d_com_dist[filt],
                                            weights2=dW[filt],verbose=verbose,is_comoving_dist=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                DD_counts[i,:] /= np.sum(dW)**2 # normalize by product of sum of weights
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
                # NB: these are already normalized 
        else:
            # Compute RR pair counts
            for i,j in enumerate(J_regions):
                # Compute pair counts between jackknife region and entire survey regions
                print("Computing RR pair counts for non-empty jackknife %d of %d"%(i+1,N_jack))
                if cross_term:
                    filt = np.where(rJ_RR==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX2_RR,rY2_RR,rZ2_RR,weights1=rW2_RR,
                                                weight_type='pair_product',X2=rX_RR[filt],Y2=rY_RR[filt],Z2=rZ_RR[filt],
                                                weights2=rW_RR[filt],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=0.5*cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    filt2 = np.where(rJ2_RR==j)
                    if len(filt2[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX_RR,rY_RR,rZ_RR,weights1=rW_RR,
                                           weight_type='pair_product',X2=rX2_RR[filt2],Y2=rY2_RR[filt2],Z2=rZ2_RR[filt2],
                                           weights2=rW2_RR[filt2],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=0.5*cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                    RR_counts[i,:] /= np.sum(rW_RR)*np.sum(rW2_RR) # normalize by product of sum of weights
                else:
                    filt = np.where(rJ_RR==j)
                    if len(filt[0])>0:
                        cross_RR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,rX_RR,rY_RR,rZ_RR,weights1=rW_RR,
                                                weight_type='pair_product',X2=rX_RR[filt],Y2=rY_RR[filt],Z2=rZ_RR[filt],
                                                weights2=rW_RR[filt],verbose=verbose,periodic=True)
                        RR_counts[i,:]+=cross_RR[:]['npairs']*cross_RR[:]['weightavg']
                        RR_counts[i,:] /= np.sum(rW_RR)**2 # normalize by product of sum of weights
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
                                             weights1=dW[filt],weight_type='pair_product',X2=rX2_DR,Y2=rY2_DR,
                                             Z2 = rZ2_DR, weights2 = rW2_DR, verbose=verbose,periodic=True)
                    D1R2_counts[i,:] +=0.5*cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                filt2 = np.where(rJ2_DR==j)
                if len(filt2[0])>0:
                    cross_D1R2 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                             weights1=dW,weight_type='pair_product',X2=rX2_DR[filt2],Y2=rY2_DR[filt2],
                                             Z2 = rZ2_DR[filt2], weights2 = rW2_DR[filt2], verbose=verbose,periodic=True)
                    D1R2_counts[i,:] +=0.5*cross_D1R2[:]['npairs']*cross_D1R2[:]['weightavg']
                D1R2_counts[i,:] /= np.sum(dW)*np.sum(rW2_DR) # normalize by product of sum of weights
                
                ## Compute D2R1 term
                filt = np.where(dJ2==j)
                if len(filt[0])>0:
                    cross_D2R1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX2[filt],dY2[filt],dZ2[filt],
                                             weights1=dW2[filt],weight_type='pair_product',X2=rX_DR,Y2=rY_DR,
                                             Z2 = rZ_DR, weights2 = rW_DR, verbose=verbose,periodic=True)
                    D2R1_counts[i,:] +=0.5*cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                filt2 = np.where(rJ_DR==j)
                if len(filt2[0])>0:
                    cross_D2R1 = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX2,dY2,dZ2,
                                             weights1=dW2,weight_type='pair_product',X2=rX_DR[filt2],Y2=rY_DR[filt2],
                                             Z2 = rZ_DR[filt2], weights2 = rW_DR[filt2], verbose=verbose,periodic=True)
                    D2R1_counts[i,:] +=0.5*cross_D2R1[:]['npairs']*cross_D2R1[:]['weightavg']
                D2R1_counts[i,:] /= np.sum(dW2)*np.sum(rW_DR) # normalize by product of sum of weights
            else:
                
                # Just compute DR term
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt],dY[filt],dZ[filt],
                                        weights1=dW[filt],weight_type='pair_product', X2=rX_DR, Y2=rY_DR, 
                                        Z2 = rZ_DR, weights2 = rW_DR, verbose=verbose,periodic=True)
                    DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                filt2 = np.where(rJ_DR==j)
                if len(filt2[0])>0:
                    cross_DR = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,
                                        weights1=dW,weight_type='pair_product', X2=rX_DR[filt2], Y2=rY_DR[filt2], 
                                        Z2 = rZ_DR[filt2], weights2 = rW_DR[filt2], verbose=verbose,periodic=True)
                    DR_counts[i,:] += 0.5*cross_DR[:]['npairs']*cross_DR[:]['weightavg']
                DR_counts[i,:] /= np.sum(dW)*np.sum(rW_DR) # normalize by product of sum of weights
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
                    DD_counts[i,:]+=0.5*cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                filt2 = np.where(dJ==j)
                if len(filt2[0])>0:
                    cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX[filt2],dY[filt2],dZ[filt2],weights1=dW[filt2],
                                            weight_type='pair_product',X2=dX2,Y2=dY2,Z2=dZ2,
                                            weights2=dW2,verbose=verbose,periodic=True)
                    DD_counts[i,:]+=0.5*cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                DD_counts[i,:] /= np.sum(dW)*np.sum(dW2) # normalize by product of sum of weights
            else:
                filt = np.where(dJ==j)
                if len(filt[0])>0:
                    cross_DD = DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,
                                            weight_type='pair_product',X2=dX[filt],Y2=dY[filt],Z2=dZ[filt],
                                            weights2=dW[filt],verbose=verbose,periodic=True)
                    DD_counts[i,:]+=cross_DD[:]['npairs']*cross_DD[:]['weightavg']
                DD_counts[i,:] /= np.sum(dW)**2 # normalize by product of sum of weights
        print("Finished DD pair counts after %d seconds"%(time.time()-init))
    
    xi_function = np.zeros_like(RR_counts)
       
    ## Now compute correlation function
    if cross_term:
        # Now use Landay-Szelay estimator:
        for j in range(N_jack):
            xi_function[j] = DD_counts[j]/RR_counts[j] - D1R2_counts[j]/RR_counts[j] - D2R1_counts[j]/RR_counts[j] + 1.
    else:
        # Now use Landay-Szelay estimator:
        for j in range(N_jack):
            xi_function[j] = DD_counts[j]/RR_counts[j] - 2.*DR_counts[j]/RR_counts[j] + 1.

    return xi_function

## Now compute correlation functions.

print("\nCOMPUTING xi_11 CORRELATION\n")
xi_11=compute_xi(random1_RR,random1_DR,data1,N_rand1_RR,N_rand1_DR,N_gal1,cross_term=False,RRname=RRname11,verbose=False)
print("\nCOMPUTING xi_12 CORRELATION\n")
xi_12=compute_xi(random1_RR,random1_DR,data1,N_rand1_RR,N_rand1_DR,N_gal1,random2_RR,random2_DR,data2,N_rand2_RR,N_rand2_DR,N_gal2,cross_term=True,RRname=RRname12,verbose=False)
print("\nCOMPUTING xi_22 CORRELATION\n")
xi_22=compute_xi(random2_RR,random2_DR,data2,N_rand2_RR,N_rand2_DR,N_gal2,cross_term=False,RRname=RRname22,verbose=False)

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
    with open(os.path.join(outdir, outname), "w+") as outfile:
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

print("NB: Number of galaxies in field 1 is %d"%N_gal1)

print("NB: Number of galaxies in field 2 is %d"%N_gal2)
