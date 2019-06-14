### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a single survey. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import sys
import numpy as np
import scipy.spatial as ss
from scipy.optimize import curve_fit

# PARAMETERS
if (len(sys.argv)!=5) and (len(sys.argv)!=6):
    print("Usage: python compute_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS}] ")
    sys.exit()
gal_file = str(sys.argv[1])
binfile = str(sys.argv[2])
outdir=str(sys.argv[3])
periodic = int(sys.argv[4])
if periodic:
    print("Assuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere");
   
else:
    RR_file = str(sys.argv[5])

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
    nw2_bar = n_bar**2*w_bar**2
else:
    nw2_bar = np.mean(gal_n**2*gal_w**2)

# Load in binning files 
r_bins = np.loadtxt(binfile)
n=len(r_bins)

## Define normalization constant
norm = V*nw2_bar

if periodic:
    ## Output periodic survey correction function
    fit_params = np.zeros([n,7])
    fit_params[:,0] = 1
    fit_params[:,3] = 1

else:
    ## load in RR counts
    RR_flat = np.loadtxt(RR_file)*np.sum(gal_w)**2. # change normalization here
    m=len(RR_flat)//n
    RR_true = RR_flat.reshape((n,m))

    # Find binning centers
    r_cen = np.mean(r_bins,axis=1)

    # Find volume of each bin
    vol_r = 4.*np.pi/3.*(r_bins[:,1]**3-r_bins[:,0]**3)

    mu_cen = np.arange(1/(2*m),1.+1/(2*m),1/m)
    delta_mu = mu_cen[-1]-mu_cen[-2]
    assert(m==len(mu_cen))

    ## Define RR model
    def RR_model(r_bin,mu):
        return vol_r[r_bin]*V*nw2_bar*delta_mu

    # Compute correction functions
    Phi_values = []
    for r_bin in range(n):
        Phi_values.append(RR_model(r_bin,mu_cen)/RR_true[r_bin,:])

    ## check for order of magnitude consistency
    if np.mean(Phi_values)>10:
        raise Exception("RR_true seems to be much smaller than RR_model. Is the input RR correctly normalized?")
    if np.mean(Phi_values)<0.1:
        raise Exception("RR_true seems to be much larger than RR_model. Is the input RR correctly normalized?")

    ## Define Phi model (piecewise continuous polynomial)
    mu_crit=0.75
    def Phi_model(mu,a0,a1,a2,b2,b3,mu_crit=mu_crit):
        b1=a1+2*mu_crit*(a2-b2)-3*b3*mu_crit**2.
        b0=a0+(a1-b1)*mu_crit+(a2-b2)*mu_crit**2.-b3*mu_crit**3.
        filt1=np.where(mu<mu_crit)
        filt2=np.where(mu>=mu_crit)
        output=np.zeros_like(mu)
        output[filt1]=a0+a1*mu[filt1]+a2*mu[filt1]**2.
        output[filt2]=b0+b1*mu[filt2]+b2*mu[filt2]**2.+b3*mu[filt2]**3.
        return output

    ## Find optimal parameters
    def fit_model(mu,good_param):
        return Phi_model(mu,good_param[0],good_param[1],good_param[2],good_param[3],good_param[4])

    fit_params=[]
    errors=[]
    for i in range(n):
        good_param,_=curve_fit(Phi_model,mu_cen,Phi_values[i],p0=[0,0,0,0,0])
        a0,a1,a2,b2,b3=good_param
        b1=a1+2*mu_crit*(a2-b2)-3*b3*mu_crit**2.
        b0=a0+(a1-b1)*mu_crit+(a2-b2)*mu_crit**2.-b3*mu_crit**3.
        out_params=[a0,a1,a2,b0,b1,b2,b3]
        fit_params.append(out_params)
        errors.append(np.abs(Phi_values[i]-fit_model(mu_cen,good_param))/fit_model(mu_cen,good_param))
            
    fit_params = np.asarray(fit_params)

    print("\nFitted with mean fractional error %.1e"%np.mean(errors))

if periodic:
    outfile = outdir+'BinCorrectionFactor_n%d_periodic_11.txt'%(n)
else:
    outfile = outdir+'BinCorrectionFactor_n%d_m%d_11.txt'%(n,m)
    
with open(outfile,"w+") as out:
    for i in range(n):
        for j in range(7):
            out.write("%.8e"%(fit_params[i,j]/norm))
            if j<6:
                out.write("\t")
            else:
                out.write("\n")     
print("\nSaved (normalized) output to %s\n"%outfile)
