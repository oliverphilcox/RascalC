## Script to compute estimate of the correlation function(s) xi(r,mu) for an entire survey via the DD/RR-1 estimator for one/two sets of tracer particles (galaxies).
## This assumes a periodic simulation, with the random RR counts computed analytically.
## Mu (RSD angle) is measured from the Z-axis.
## If two sets of galaxies are provided, the code will automatically compute xi_11 xi_12 and xi_22, else only xi_11 will be computed.

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=8:
    if len(sys.argv)!=9:
        print("Usage: python xi_estimator_periodic.py {GALAXY_FILE} {RADIAL_BIN_FILE} {BOXSIZE} {MU_MAX} {N_MU_BINS} {NTHREADS} {OUTPUT_DIR} [{GALAXY_FILE_2}]")
        sys.exit(1)
Dname = str(sys.argv[1])
binfile = str(sys.argv[2])
boxsize = float(sys.argv[3])
mu_max = float(sys.argv[4])
nmu_bins = int(sys.argv[5])
nthreads = int(sys.argv[6])
outdir=str(sys.argv[7])

if len(sys.argv)==9:
    multifield = True
    Dname2 = str(sys.argv[8])
else:
    multifield = False

## First read in weights and positions:
dtype = np.double

print("Counting lines in galaxy file")
total_lines=0
for n, line in enumerate(open(Dname, 'r')):
    total_lines+=1

dX,dY,dZ,dW=[np.zeros(total_lines) for _ in range(4)]

print("Reading in galaxy data");
for n, line in enumerate(open(Dname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(), dtype=float)
    dX[n]=split_line[0];
    dY[n]=split_line[1];
    dZ[n]=split_line[2];
    if len(split_line)>3:
        dW[n]=split_line[3];
    else:
        dW[n]=1.

N_gal = len(dX) # number of particles

print("Number of galaxy particles: %.1e"%N_gal)

## Check for periodicity
xrange = max(dX)-min(dX)
yrange = max(dY)-min(dY)
zrange = max(dZ)-min(dZ)

assert(np.abs(xrange-boxsize)/boxsize<0.001),'Data is not periodic! X-range is %.2e compared to boxsize %.2e'%(xrange,boxsize)
assert(np.abs(yrange-boxsize)/boxsize<0.001),'Data is not periodic! Y-range is %.2e compared to boxsize %.2e'%(yrange,boxsize)
assert(np.abs(zrange-boxsize)/boxsize<0.001),'Data is not periodic! Z-range is %.2e compared to boxsize %.2e'%(zrange,boxsize)

if multifield:

    print("Counting lines in galaxy file 2")
    total_lines2=0
    for n, line in enumerate(open(Dname2,'r')):
        total_lines2+=1

    dX2,dY2,dZ2,dW2=[np.zeros(total_lines2) for _ in range(4)]

    print("Reading in galaxy data for galaxy 2");
    for n, line in enumerate(open(Dname2, 'r')):
        if n%1000000==0:
            print("Reading line %d of %d" %(n,total_lines2))
        split_line=np.array(line.split(), dtype=float)
        dX2[n]=split_line[0];
        dY2[n]=split_line[1];
        dZ2[n]=split_line[2];
        dW2[n]=split_line[3];

    N_gal2 = len(dX2) # number of particles
    print("Number of galaxy particles in second set: %.1e"%N_gal2)

    ## Check for periodicity
    xrange = max(dX2)-min(dX2)
    yrange = max(dY2)-min(dY2)
    zrange = max(dZ2)-min(dZ2)

    assert(np.abs(xrange-boxsize)/boxsize<0.001),'Data set 2 is not periodic! X-range is %.2e compared to boxsize %.2e'%(xrange,boxsize)
    assert(np.abs(yrange-boxsize)/boxsize<0.001),'Data set 2 is not periodic! Y-range is %.2e compared to boxsize %.2e'%(yrange,boxsize)
    assert(np.abs(zrange-boxsize)/boxsize<0.001),'Data set 2 is not periodic! Z-range is %.2e compared to boxsize %.2e'%(zrange,boxsize)

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

print("Using periodic input data");
from Corrfunc.theory.DDsmu import DDsmu

# Compute RR counts analytically
RR_counts = 4.*np.pi/3.*(all_bins[:,1]**3-all_bins[:,0]**3)/boxsize**3*mu_max/nmu_bins

import time
init = time.time()

# Now compute DD counts
print("Computing DD pair counts")
tmpDD=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',verbose=True,periodic=True)
DD_counts = tmpDD[:]['npairs']*tmpDD[:]['weightavg']
DD_counts/=np.sum(dW)**2.
print("Finished after %d seconds"%(time.time()-init))

# Now use CF estimator:
xi_reshape = DD_counts.reshape(nrbins,nmu_bins)/RR_counts.reshape(nrbins,1) - 1.

if multifield:
    # Compute cross fields
    init  = time.time()
    print("Computing DD pair counts for cross pair counts")
    tmpDD12=DDsmu(0,nthreads,binfile,mu_max,nmu_bins,dX,dY,dZ,weights1=dW,weight_type='pair_product',
                X2=dX2,Y2=dY2,Z2=dZ2,weights2=dW2,
                verbose=True,periodic=True)
    DD12_counts = tmpDD12[:]['npairs']*tmpDD12[:]['weightavg']
    DD12_counts/=np.sum(dW)*np.sum(dW2)
    print("Finished after %d seconds"%(time.time()-init))

    xi_reshape12 = DD12_counts.reshape(nrbins,nmu_bins)/RR_counts.reshape(nrbins,1)-1.

    # Compute second field pair counts
    init = time.time()
    print("Computing DD pair counts for second dataset")
    tmpDD2=DDsmu(1,nthreads,binfile,mu_max,nmu_bins,dX2,dY2,dZ2,weights1=dW2,weight_type='pair_product',
                verbose=True,periodic=True)
    DD2_counts = tmpDD2[:]['npairs']*tmpDD2[:]['weightavg']
    DD2_counts/=np.sum(dW2)**2.
    print("Finished after %d seconds"%(time.time()-init))

    xi_reshape2 = DD2_counts.reshape(nrbins,nmu_bins)/RR_counts.reshape(nrbins,1)-1.

# Save output files:
import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Define mu centers
mean_mus = np.linspace(0.5/nmu_bins,1-0.5/nmu_bins,nmu_bins)

outname='xi_n%d_m%d_periodic_11.dat'%(nrbins,nmu_bins)
print("Saving correlation function(s)")
with open(os.path.join(outdir, outname), "w+") as outfile:
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

if multifield:

    outname='xi_n%d_m%d_periodic_12.dat'%(nrbins,nmu_bins)
    with open(os.path.join(outdir, outname), "w+") as outfile:
        for r in mean_bins:
            outfile.write("%.8e "%r)
        outfile.write("\n")
        for mu in mean_mus:
            outfile.write("%.8e "%mu)
        outfile.write("\n");
        for i in range(nrbins):
            for j in range(nmu_bins):
                outfile.write("%.8e "%xi_reshape12[i,j])
            outfile.write("\n")

    print("Cross correlation function written successfully to %s"%(outdir+outname))


    outname='xi_n%d_m%d_periodic_22.dat'%(nrbins,nmu_bins)
    with open(os.path.join(outdir, outname), "w+") as outfile:
        for r in mean_bins:
            outfile.write("%.8e "%r)
        outfile.write("\n")
        for mu in mean_mus:
            outfile.write("%.8e "%mu)
        outfile.write("\n");
        for i in range(nrbins):
            for j in range(nmu_bins):
                outfile.write("%.8e "%xi_reshape2[i,j])
            outfile.write("\n")

    print("Second field correlation function written successfully to %s"%(outdir+outname))

if multifield:
    print("NB: Number of galaxies is (%d, %d)"%(N_gal,N_gal2))
if multifield:
    print("NB: Number of galaxies is %d"%N_gal)
