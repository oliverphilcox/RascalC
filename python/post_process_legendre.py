## Script to post-process the single-field Legendre binned integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (6, 7, 8):
    print("Usage: python post_process_legendre.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING} [{SKIP_R_BINS}]]")
    sys.exit()

file_root = str(sys.argv[1])
n = int(sys.argv[2])
max_l = int(sys.argv[3])
n_samples = int(sys.argv[4])
outdir = str(sys.argv[5])
alpha = float(sys.argv[6]) if len(sys.argv) >= 7 else 1.
skip_bins = int(sys.argv[7]) if len(sys.argv) >= 8 else 0

# mask for skipping bins
r_mask = np.arange(n) >= skip_bins

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

def load_matrices(index):
    """Load intermediate or full covariance matrices"""
    cov_root = os.path.join(file_root, 'CovMatricesAll/')
    c2 = np.loadtxt(cov_root+'c2_n%d_l%d_11_%s.txt'%(n,max_l,index))
    c3 = np.loadtxt(cov_root+'c3_n%d_l%d_1,11_%s.txt'%(n,max_l,index))
    c4 = np.loadtxt(cov_root+'c4_n%d_l%d_11,11_%s.txt'%(n,max_l,index))

    N = len(c2)
    assert N % n == 0, "Number of bins mismatch"
    nrepeat = N // n
    full_mask = np.array(nrepeat * list(r_mask)) # repeat the r_mask since cov terms are first ordered by l
    c2, c3, c4 = (a[full_mask][:, full_mask] for a in (c2, c3, c4)) # select rows and columns

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
    c2t,c3t,c4t=load_matrices(i)
    c2s.append(c2t)
    c3s.append(c3t)
    c4s.append(c4t)
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

output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_n%d_l%d.npz'%(n,max_l))
np.savez(output_name,full_theory_covariance=full_cov,
         shot_noise_rescaling=alpha,full_theory_precision=full_prec,
         N_eff=N_eff_D,full_theory_D_matrix=full_D_est,
         individual_theory_covariances=partial_cov)

print("Saved output covariance matrices as %s"%output_name)
