## Script to post-process the single-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv)!=6 and len(sys.argv)!=7:
    print("Usage: python post_process_default.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]")
    sys.exit()

file_root = str(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
n_samples = int(sys.argv[4])
outdir = str(sys.argv[5])
if len(sys.argv)==7:
    alpha = float(sys.argv[6])
else:
    alpha = 1.;

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

def load_matrices(index):
    """Load intermediate or full covariance matrices"""
    cov_root = os.path.join(file_root, 'CovMatricesAll/')
    c2 = np.diag(np.loadtxt(cov_root+'c2_n%d_m%d_11_%s.txt'%(n,m,index)))
    c3 = np.loadtxt(cov_root+'c3_n%d_m%d_1,11_%s.txt'%(n,m,index))
    c4 = np.loadtxt(cov_root+'c4_n%d_m%d_11,11_%s.txt'%(n,m,index))

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

# Compute full covariance matrices and precision
full_cov = c4+c3*alpha+c2*alpha**2.
n_bins = len(c4)

# Compute full precision matrix
print("Computing the full precision matrix estimate:")
# Load in partial theoretical matrices
c2s,c3s,c4s=[],[],[]
for i in range(n_samples):
    print("Loading full subsample %d of %d"%(i+1,n_samples))
    c2,c3,c4=load_matrices(i)
    c2s.append(c2)
    c3s.append(c3)
    c4s.append(c4)
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

output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Default_n%d_m%d.npz'%(n,m))
np.savez(output_name,full_theory_covariance=full_cov,
         shot_noise_rescaling=alpha,full_theory_precision=full_prec,
         N_eff=N_eff_D,full_theory_D_matrix=full_D_est,
         individual_theory_covariances=partial_cov)

print("Saved output covariance matrices as %s"%output_name)
