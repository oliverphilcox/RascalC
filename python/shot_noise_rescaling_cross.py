## Script to estimate the shot noise rescaling parameters using the jackknife covariance matrices of a single survey for multiple fields. 
## This saves the optimal shot noise parameters and saves the full (non-jackknife) theoretical covariance matrices. In addition, we save data and theoretical jackknife covariance matrices of each type

import numpy as np
import sys

# PARAMETERS
if len(sys.argv)!=7:
    print("Usage: python shot_noise_rescaling.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}")
    sys.exit()
        
jackknife_file = str(sys.argv[1])
weight_dir = str(sys.argv[2])
file_root = str(sys.argv[3])
m = int(sys.argv[4])
n_samples = int(sys.argv[5])
outdir = str(sys.argv[6])

# Load jackknife xi estimates from data
print("Loading correlation function jackknife estimates from %s"%jackknife_file)
xi_jack = np.loadtxt(jackknife_file,skiprows=2)
n_bins = xi_jack.shape[1] # total bins
n_jack = xi_jack.shape[0] # total jackknives
n = n_bins//m # radial bins

weight_file = weight_dir+'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack)
RR_file = weight_dir+'binned_pair_counts_n%d_m%d_j%d_11.dat'%(n,m,n_jack)

print("Loading weights file from %s"%weight_file)
weights = np.loadtxt(weight_file)[:,1:]

# First exclude any dodgy jackknife regions
good_jk=np.unique(np.where(np.isfinite(xi_jack))[0])
print("Using %d out of %d jackknives"%(n_jack,len(good_jk)))
xi_jack = xi_jack[good_jk]
weights = weights[good_jk]

# Compute data covariance matrix
print("Computing data covariance matrix")
mean_xi = np.mean(xi_jack,axis=0)
tmp = weights*(xi_jack-np.mean(xi_jack,axis=0))
data_cov = np.matmul(tmp.T,tmp)
denom = np.matmul(weights.T,weights)
data_cov /= (np.ones_like(denom)-denom)

print("Loading weights file from %s"%RR_file)
RR=np.loadtxt(RR_file)

def load_matrices(index,jack=True):
    """Load intermediate or full covariance matrices"""
    if jack:
        cov_root = file_root+'CovMatricesJack/'
    else:
        cov_root = file_root+'CovMatricesAll/'
    c2 = np.diag(np.loadtxt(cov_root+'c2_n%d_m%d_11_%s.txt'%(n,m,index)))
    c3 = np.loadtxt(cov_root+'c3_n%d_m%d_1,11_%s.txt'%(n,m,index))
    c4 = np.loadtxt(cov_root+'c4_n%d_m%d_11,11_%s.txt'%(n,m,index))
    if jack:
        EEaA1 = np.loadtxt(cov_root+'EE1_n%d_m%d_11_%s.txt' %(n,m,index))
        EEaA2 = np.loadtxt(cov_root+'EE2_n%d_m%d_11_%s.txt' %(n,m,index))
        RRaA1 = np.loadtxt(cov_root+'RR1_n%d_m%d_11_%s.txt' %(n,m,index))
        RRaA2 = np.loadtxt(cov_root+'RR2_n%d_m%d_11_%s.txt' %(n,m,index))
    
        # Compute disconnected term
        w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
        w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
        diff1 = EEaA1-w_aA1*EEaA1.sum(axis=0)
        diff2 = EEaA2-w_aA2*EEaA2.sum(axis=0)
        RRaRRb=np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
        fact=np.eye(n_bins)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
        cx = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
        c4+=cx
    
    return c2,c3,c4

# Load in full jackknife theoretical matrices
print("Loading best estimate of jackknife covariance matrix")
c2,c3,c4=load_matrices('full')

# Load in partial jackknife theoretical matrices
c2s,c3s,c4s=[],[],[]
for i in range(n_samples):
    print("Loading subsample %d of %d"%(i+1,n_samples))
    c2,c3,c4=load_matrices(i)
    c2s.append(c2)
    c3s.append(c3)
    c4s.append(c4)

# Compute inverted matrix
def Psi(alpha):
    c_tot = c2*alpha**2.+c3*alpha+c4
    partial_cov=[]
    for i in range(n_samples):
        partial_cov.append(alpha**2.*c2s[i]+alpha*c3s[i]+c4s[i])
    tmp=0.
    for i in range(n_samples):
        c_excl_i = np.mean(partial_cov[:i]+partial_cov[i+1:],axis=0)
        tmp+=np.matmul(np.linalg.inv(c_excl_i),partial_cov[i])
    D_est=(n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
    Psi = np.matmul(np.eye(n_bins)-D_est,np.linalg.inv(c_tot))
    return Psi

def neg_log_L1(alpha):
    """Return negative log L1 likelihood between data and theory covariance matrices"""
    Psi_alpha = Psi(alpha)
    logdet = np.linalg.slogdet(Psi_alpha)
    if logdet[0]<0:
        # Remove any dodgy inversions
        return np.inf        
    return np.trace(np.matmul(Psi_alpha,data_cov))-logdet[1]

# Now optimize for shot-noise rescaling parameter alpha
print("Optimizing for the shot-noise rescaling parameter")
from scipy.optimize import fmin
alpha_best = fmin(neg_log_L1,1.)
print("Optimization complete")

# Now save to file
np.save(outdir+'alpha_estimate.npy',alpha_best)

jack_cov = c4+c3*alpha_best+c2*alpha_best**2.
c2f,c3f,c4f=load_matrices('full',jack=False)
full_cov = c4f+c3f*alpha_best+c2f*alpha_best

output_name = outdir+'Rescaled_Covariance_Matrices_n%d_m%d_j%d.npz'%(n,m,n_jack)
np.savez(output_name,jackknife_theory=jack_cov,full_theory=full_cov,jackknife_data=data_cov,shot_noise_rescaling=alpha_best)
print("Saved output covariance matrices as %s",output_name)
