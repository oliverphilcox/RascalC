### Function to fit a model to the 3PCF survey correction function, defined as the ratio between model and true RR pair counts for a single survey. This fits a piecewise polynomial model to the data.

## NB: Input RRR counts should be normalized by summed cubed random weights here.
## NB: Assume mu is in [-1,1] limit here


import sys
import os
import numpy as np
import scipy.spatial as ss

# PARAMETERS
if (len(sys.argv)!=5) and (len(sys.argv)!=6):
    print("Usage python compute_3pcf_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RRR_COUNTS}]")
    sys.exit(1);

gal_file = str(sys.argv[1])
binfile = str(sys.argv[2])
outdir=str(sys.argv[3])
periodic = int(sys.argv[4])
if periodic:
    if(len(sys.argv)!=5):
        print("Usage python compute_3pcf_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RRR_COUNTS}]")
        sys.exit(1);
    print("\nAssuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere");
   
else:
    if(len(sys.argv)!=6):
        print("Usage python compute_3pcf_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RRR_COUNTS}]")
        sys.exit(1);
    RRR_file = str(sys.argv[5])

## Load galaxies
print("\nLoading galaxies")
all_gal = np.loadtxt(gal_file)
gal_x = all_gal[:,0]
gal_y = all_gal[:,1]
gal_z = all_gal[:,2]
gal_w = all_gal[:,3]
gal_n = (1./gal_w-1.)/20000.

N_gal = len(all_gal)
w_bar = np.mean(gal_w)

## Find survey volume via ConvexHull in Scipy
hull = ss.ConvexHull(np.vstack([gal_x,gal_y,gal_z]).T)
print('\nSurvey volume is approximately: %.2f (Gpc/h)^3'%(hull.volume/1e9))
V=hull.volume # in (Mpc/h)^3

## Galaxy number density
n_bar = N_gal/V

if periodic:
    nw3_bar = n_bar**3*w_bar**3
else:
    nw3_bar = np.mean(gal_n**3*gal_w**3)

# Load in binning files 
r_bins = np.loadtxt(binfile)
n=len(r_bins)

## Define normalization constant
norm = 6.*V*nw3_bar

print("Normalizing output survey correction by %.2e"%norm);

if periodic:
    
    ## Output periodic survey correction function
    phi_inv_mult = np.zeros([n,n,7]);
    
    ## Set to correct periodic survey values
    phi_inv_mult[:,:,0]=1.

else:
    from scipy.special import legendre
    
    ## Load triple counts and renormalize
    tmp_triple_counts = np.loadtxt(RRR_file)*np.sum(gal_w)**3
    
    # Compute number of angular bins in data-set
    m = (len(tmp_triple_counts)//n)//n
    assert len(tmp_triple_counts)%m==0, "Incorrect RRR format"
    

    mu_all = np.linspace(-1,1,m+1)
    mu_cen = 0.5*(mu_all[1:]+mu_all[:-1])
    
    RRR_true = np.zeros([n,n,m])
    
    ## load in RRR counts (and add symmetries)
    for i in range(len(tmp_triple_counts)):
        RRR_true[(i//m)//n,(i//m)%n,i%m] += tmp_triple_counts[i]*0.5
        RRR_true[(i//m)%n,(i//m)//n,i%m] += tmp_triple_counts[i]*0.5
        
    ## Now construct Legendre moments
    leg_triple = np.zeros([n,n,7])
    for a in range(n):
        for b in range(n):
            for ell in range(7):
                # (NB: we've absorbed a factor of delta_mu into RRR_true here)
                leg_triple[a,b,ell]+=np.sum(legendre(ell)(mu_cen)*RRR_true[a,b,:])*(2.*ell+1.)

    vol_r = lambda b: 4.*np.pi/3.*(r_bins[b,1]**3.-r_bins[b,0]**3.)

    ## Construct inverse multipoles of Phi
    phi_inv_mult = np.zeros([n,n,7])
    for b1 in range(n):
        for b2 in range(n):
            phi_inv_mult[b1,b2,:] = leg_triple[b1,b2,:]/(3.*nw3_bar*V*vol_r(b1)*vol_r(b2))
            
    ## Check all seems reasonable
    if np.mean(phi_inv_mult[:,:,0])<1e-3:
        print(phi_inv_mult[:,:,0])
        print("Survey correction function seems too small - are the RRR counts normalized correctly?")
        sys.exit(1)
    if np.mean(phi_inv_mult[:,:,0])>1e3:
        print("Survey correction function seems too large - are the RRR counts normalized correctly?")
        sys.exit(1)
    
if periodic:
    outfile = os.path.join(outdir, 'BinCorrectionFactor3PCF_n%d_periodic.txt'%(n))
else:
    outfile = os.path.join(outdir, 'BinCorrectionFactor3PCF_n%d_m%d.txt'%(n,m))
    
with open(outfile,"w+") as out:
    for b1 in range(n):
        for b2 in range(n):
            for ell in range(7):
                out.write("%.8e"%(phi_inv_mult[b1,b2,ell]*norm))
                if ell<6:
                    out.write("\t")
                if ell==6:
                    out.write("\n")
print("\nSaved (normalized) output to %s\n"%outfile)
