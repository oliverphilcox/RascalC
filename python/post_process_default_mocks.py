## Script to post-process the single-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (7, 8):
    print("Usage: python post_process_default_mocks.py {MOCK_COV_FILE} {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SKIP_R_BINS}]")
    sys.exit()
        
mock_cov_file = str(sys.argv[1])
file_root = str(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])
n_samples = int(sys.argv[5])
outdir = str(sys.argv[6])
skip_bins = 0
if len(sys.argv) == 8:
    skip_bins = int(sys.argv[7]) * m # convert from radial to total number of bins right away

mock_cov = np.loadtxt(mock_cov_file)[skip_bins:, skip_bins:] # load external mock covariance matrix

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

def load_matrices(index):
    """Load intermediate or full covariance matrices"""
    cov_root = os.path.join(file_root, 'CovMatricesAll/')
    c2 = np.diag(np.loadtxt(cov_root+'c2_n%d_m%d_11_%s.txt'%(n,m,index))[skip_bins:])
    c3 = np.loadtxt(cov_root+'c3_n%d_m%d_1,11_%s.txt'%(n,m,index))[skip_bins:, skip_bins:]
    c4 = np.loadtxt(cov_root+'c4_n%d_m%d_11,11_%s.txt'%(n,m,index))[skip_bins:, skip_bins:]

    # Now symmetrize and return matrices
    return c2,0.5*(c3+c3.T),0.5*(c4+c4.T)

# Load in full theoretical matrices
print("Loading best estimate of covariance matrix")
c2,c3,c4=load_matrices('full')

# Check matrix convergence
from numpy.linalg import eigvalsh
eig_c4 = eigvalsh(c4)
eig_c2 = eigvalsh(c2)
if min(eig_c4)<-1.*min(eig_c2):
    print("4-point covariance matrix has not converged properly via the eigenvalue test. Exiting")
    print("Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4), min(eig_c2)))
    sys.exit()

n_bins = len(c4)

# Load in partial theoretical matrices
c2s,c3s,c4s=[],[],[]
for i in range(n_samples):
    print("Loading full subsample %d of %d"%(i+1,n_samples))
    c2,c3,c4=load_matrices(i)
    c2s.append(c2)
    c3s.append(c3)
    c4s.append(c4)

# Compute inverted matrix
def Psi(alpha):
    """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
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
    return np.trace(np.matmul(Psi_alpha, mock_cov))-logdet[1]

# Now optimize for shot-noise rescaling parameter alpha
print("Optimizing for the shot-noise rescaling parameter")
from scipy.optimize import fmin
alpha_best = fmin(neg_log_L1,1.)
print("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)
alpha = alpha_best # to save editing later

# Compute full covariance matrices and precision
full_cov = c4+c3*alpha+c2*alpha**2.

# Compute full precision matrix
print("Computing the full precision matrix estimate:")
partial_cov=[]
for i in range(n_samples):
    partial_cov.append(alpha**2.*c2s[i]+alpha*c3s[i]+c4s[i])
tmp=0.
for i in range(n_samples):
    c_excl_i = np.mean(partial_cov[:i]+partial_cov[i+1:],axis=0)
    tmp+=np.matmul(np.linalg.inv(c_excl_i),partial_cov[i])
full_D_est=(n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
full_prec = np.matmul(np.eye(n_bins)-full_D_est,np.linalg.inv(full_cov))
print("Full precision matrix estimate computed")

# Now compute effective N:
slogdetD=np.linalg.slogdet(full_D_est)
D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
if slogdetD[0]<0:
    print("N_eff is negative! Setting to zero")
    N_eff_D = 0.
else:
    N_eff_D = (n_bins+1.)/D_value+1.
    print("Total N_eff Estimate: %.4e"%N_eff_D)

output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Default_Mocks_n%d_m%d.npz'%(n,m))
np.savez(output_name, full_theory_covariance=full_cov,
         shot_noise_rescaling=alpha, full_theory_precision=full_prec,
         N_eff=N_eff_D, full_theory_D_matrix=full_D_est,
         individual_theory_covariances=partial_cov, mock_covariance=mock_cov)

print("Saved output covariance matrices as %s"%output_name)
