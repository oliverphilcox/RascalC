### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a two surveys. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import sys
import numpy as np
import scipy.spatial as ss
from scipy.optimize import curve_fit

# PARAMETERS
if len(sys.argv)!=8:
    print("Usage: python compute_correction_function.py {GALAXY_FILE} {GALAXY_FILE_2} {BIN_FILE} {RR_COUNTS_11} {RR_COUNTS_12} {RR_COUNTS_22} {OUTPUT_DIR}")
    sys.exit()
gal_file = str(sys.argv[1])
gal_file2 = str(sys.argv[2])
binfile = str(sys.argv[3])
RR_file = str(sys.argv[4])
RR_file12 = str(sys.argv[5])
RR_file2 = str(sys.argv[6])
outdir=str(sys.argv[7])

## Load galaxies
print("\nLoading galaxy set 1")
all_gal = np.loadtxt(gal_file)

print("Loading galaxy set 2")
all_gal2 = np.loadtxt(gal_file2)

# Define survey galaxy coordinates
gal_x = all_gal[:,0]
gal_y = all_gal[:,1]
gal_z = all_gal[:,2]
gal_w = all_gal[:,3]
gal2_x = all_gal2[:,0]
gal2_y = all_gal2[:,1]
gal2_z = all_gal2[:,2]
gal2_w = all_gal2[:,3]

w_bar1 = np.mean(gal_w)
w_bar2 = np.mean(gal2_w)

N_gal = len(all_gal)
N_gal2 = len(all_gal2)

## Find survey volume via ConvexHull in Scipy
hull = ss.ConvexHull(np.vstack([gal_x,gal_y,gal_z]).T)
print('\nSurvey volume 1 is approximately: %.2f (Gpc/h)^3'%(hull.volume/1e9))
V=hull.volume # in (Mpc/h)^3

hull = ss.ConvexHull(np.vstack([gal2_x,gal2_y,gal2_z]).T)
print('Survey volume 2 is approximately: %.2f (Gpc/h)^3'%(hull.volume/1e9))
V2=hull.volume # in (Mpc/h)^3

## Galaxy number density
n_bar1 = N_gal/V
n_bar2 = N_gal2/V2

# Load in binning files 
r_bins = np.loadtxt(binfile)
n=len(r_bins)

# Find binning centers
r_cen = np.mean(r_bins,axis=1)
delta_r = r_cen[-1]-r_cen[-2]
    
# Load RR counts
RR_flat = np.loadtxt(RR_file)*np.sum(gal_w)**2. # change normalization here
m=len(RR_flat)//n
RR_true = RR_flat.reshape((n,m))

RR_flat2 = np.loadtxt(RR_file2)*np.sum(gal2_w)**2.
assert((len(RR_flat2)//n)==m), "Need same bins for all RR files."
RR_true2 = RR_flat2.reshape((n,m))

RR_flat12 = np.loadtxt(RR_file2)*np.sum(gal_w)*np.sum(gal2_w)
assert((len(RR_flat12)//n)==m), "Need same bins for all RR files."
RR_true12 = RR_flat12.reshape((n,m))

# Load mu bins
mu_cen = np.arange(1/(2*m),1.+1/(2*m),1/m)
delta_mu = mu_cen[-1]-mu_cen[-2]
assert(m==len(mu_cen))    
    
def compute_phi(w_bar1,w_bar2,n_bar1,n_bar2,V,this_RR,index):
    print("\nComputing survey correction factor %s"%index)
    ## Define normalization constant
    norm = 4.*np.pi*V*(n_bar1*w_bar1)*(n_bar2*w_bar2)
    
    ## Define RR model
    def RR_model(r_cen,mu):
        return norm*r_cen**2.*delta_r*delta_mu

    # Compute correction functions
    Phi_values = []
    for r_bin in range(n):
        Phi_values.append(RR_model(r_cen[r_bin],mu_cen)/this_RR[r_bin,:])

    ## check for order of magnitude consistency
    if np.mean(Phi_values)>50:
        raise Exception("RR_true seems to be much smaller than RR_model. Is the input RR correctly normalized?")
    if np.mean(Phi_values)<0.02:
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

    print("Fitted dataset %s with mean fractional error %.1e"%(index,np.mean(errors)))

    outfile = outdir+'BinCorrectionFactor_n%d_m%d_%s.txt'%(n,m,index)
    with open(outfile,"w+") as out:
        for i in range(n):
            for j in range(7):
                out.write("%.8e"%(fit_params[i,j]/norm))
                if j<6:
                    out.write("\t")
                else:
                    out.write("\n")    
    print("Saved (normalized) output to %s"%outfile)

# Now run this
compute_phi(w_bar1,w_bar1,n_bar1,n_bar1,V,RR_true,"11");
compute_phi(w_bar1,w_bar2,n_bar1,n_bar2,V,RR_true,"12"); # just use volume 1 here
compute_phi(w_bar2,w_bar2,n_bar2,n_bar2,V2,RR_true,"22");
