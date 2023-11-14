## Script to post-process the single-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import numpy as np
import sys, os

# PARAMETERS
if len(sys.argv) != 8:
    print("Usage: python post_process_legendre_mix_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR}")
    sys.exit(1)
        
jackknife_file = str(sys.argv[1])
weight_dir = str(sys.argv[2])
file_root = str(sys.argv[3])
m = int(sys.argv[4])
max_l = int(sys.argv[5])
assert max_l % 2 == 0, "Only even multipoles supported"
n_l = max_l//2 + 1 # number of multipoles
n_samples = int(sys.argv[6])
outdir = str(sys.argv[7])

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Load jackknife xi estimates from data
print("Loading correlation function jackknife estimates from %s" % jackknife_file)
xi_jack = np.loadtxt(jackknife_file, skiprows=2)
n_bins_smu = xi_jack.shape[1] # total s,mu bins
n_jack = xi_jack.shape[0] # total jackknives
n = n_bins_smu // m # radial bins
n_bins = n_l * n # total Legendre bins

weight_file = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_11.dat' % (n, m, n_jack))
mu_bin_legendre_file = os.path.join(weight_dir, 'mu_bin_legendre_factors_m%d_l%d.txt' % (m, max_l))

print("Loading jackknife weights from %s" % weight_file)
weights = np.loadtxt(weight_file)[:, 1:]

# First exclude any dodgy jackknife regions
good_jk = np.where(np.all(np.isfinite(xi_jack), axis=1))[0] # all xi in jackknife have to be normal numbers
print("Using %d out of %d jackknives" % (len(good_jk), n_jack))
xi_jack = xi_jack[good_jk]
weights = weights[good_jk]
weights /= np.sum(weights, axis=0) # renormalize weights after possibly discarding some jackknives

# Compute data covariance matrix
print("Computing data covariance matrix")
mean_xi = np.sum(xi_jack*weights, axis=0)
tmp = weights*(xi_jack-mean_xi)
data_cov = np.matmul(tmp.T, tmp)
denom = np.matmul(weights.T, weights)
data_cov /= (np.ones_like(denom)-denom)

print("Loading mu bin Legendre factors from %s" % mu_bin_legendre_file)
mu_bin_legendre_factors = np.loadtxt(mu_bin_legendre_file) # rows correspond to mu bins, columns to multipoles

# Project the data jackknife covariance from mu bins to Legendre multipoles
data_cov = data_cov.reshape(n, m, n, m) # make the array 4D with [r_bin, mu_bin] indices for rows and columns
data_cov = np.einsum("imjn,mp,nq->ipjq", data_cov, mu_bin_legendre_factors, mu_bin_legendre_factors) # use mu bin Legendre factors to project mu bins into Legendre multipoles, staying within the same radial bins. The indices are now [r_bin, ell] for rows and columns
data_cov = data_cov.reshape(n_bins, n_bins)

def load_matrices(index,jack=True):
    """Load intermediate or full covariance matrices"""
    if jack:
        cov_root = os.path.join(file_root, 'CovMatricesJack/')
    else:
        cov_root = os.path.join(file_root, 'CovMatricesAll/')
    c2 = np.loadtxt(cov_root+'c2_n%d_l%d_11_%s.txt' % (n, max_l, index))
    c3 = np.loadtxt(cov_root+'c3_n%d_l%d_1,11_%s.txt' % (n, max_l, index))
    c4 = np.loadtxt(cov_root+'c4_n%d_l%d_11,11_%s.txt' % (n, max_l, index))
    # Now symmetrize and return matrices
    return c2, 0.5*(c3+c3.T), 0.5*(c4+c4.T)

# Load in full jackknife theoretical matrices
print("Loading best estimate of jackknife covariance matrix")
c2j, c3j, c4j = load_matrices('full')

# Check matrix convergence
from numpy.linalg import eigvalsh
eig_c4 = eigvalsh(c4j)
eig_c2 = eigvalsh(c2j)
if min(eig_c4) < -1. * min(eig_c2):
    print("Jackknife 4-point covariance matrix has not converged properly via the eigenvalue test. Exiting")
    print("Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4), min(eig_c2)))
    sys.exit(1)

# Load in partial jackknife theoretical matrices
c2s, c3s, c4s = [], [], []
for i in range(n_samples):
    print("Loading jackknife subsample %d of %d" % (i+1, n_samples))
    c2, c3, c4=load_matrices(i)
    c2s.append(c2)
    c3s.append(c3)
    c4s.append(c4)
c2s, c3s, c4s = [np.array(a) for a in (c2s, c3s, c4s)]

# Compute inverted matrix
def Psi(alpha):
    """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
    c_tot = c2j*alpha**2.+c3j*alpha+c4j
    partial_cov = alpha**2 * c2s + alpha * c3s + c4s
    sum_partial_cov = np.sum(partial_cov, axis=0)
    tmp = 0.
    for i in range(n_samples):
        c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
        tmp += np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
    D_est = (n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
    Psi = np.matmul(np.eye(n_bins)-D_est, np.linalg.inv(c_tot))
    return Psi

def neg_log_L1(alpha):
    """Return negative log L1 likelihood between data and theory covariance matrices"""
    Psi_alpha = Psi(alpha)
    logdet = np.linalg.slogdet(Psi_alpha)
    if logdet[0] < 0:
        # Remove any dodgy inversions
        return np.inf        
    return np.trace(np.matmul(Psi_alpha,data_cov))-logdet[1]

# Now optimize for shot-noise rescaling parameter alpha
print("Optimizing for the shot-noise rescaling parameter")
from scipy.optimize import fmin
alpha_best = fmin(neg_log_L1, 1.)
print("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)

# Compute jackknife and full covariance matrices
jack_cov = c4j + c3j*alpha_best + c2j*alpha_best**2.
jack_prec = Psi(alpha_best)
c2f, c3f, c4f=load_matrices('full', jack=False)
full_cov = c4f + c3f*alpha_best + c2f*alpha_best**2.

# Check convergence
eig_c4f = eigvalsh(c4f)
eig_c2f = eigvalsh(c2f)
if min(eig_c4f)<min(eig_c2f)*-1.:
    print("Full 4-point covariance matrix has not converged properly via the eigenvalue test. Exiting")
    print("Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4f), min(eig_c2f)))
    sys.exit(1)

# Compute full precision matrix
print("Computing the full precision matrix estimate:")
# Load in partial jackknife theoretical matrices
c2fs, c3fs, c4fs = [], [], []
for i in range(n_samples):
    print("Loading full subsample %d of %d" % (i+1, n_samples))
    c2,c3,c4=load_matrices(i, jack=False)
    c2fs.append(c2)
    c3fs.append(c3)
    c4fs.append(c4)
c2fs, c3fs, c4fs = [np.array(a) for a in (c2fs, c3fs, c4fs)]
partial_cov = alpha_best**2 * c2fs + alpha_best * c3fs + c4fs
sum_partial_cov = np.sum(partial_cov, axis=0)
tmp = 0.
for i in range(n_samples):
    c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
    tmp += np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
full_D_est = (n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
full_prec = np.matmul(np.eye(n_bins)-full_D_est,np.linalg.inv(full_cov))
print("Full precision matrix estimate computed")    

# Now compute effective N:
slogdetD = np.linalg.slogdet(full_D_est)
D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
N_eff_D = (n_bins + 1.)/D_value + 1.
print("Total N_eff Estimate: %.4e" % N_eff_D)

# Jackknife covariance for posterity
partial_jack_cov = alpha_best**2 * c2s + alpha_best * c3s + c4s

output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_Jackknife_n%d_l%d_j%d.npz' % (n, max_l, n_jack))
np.savez(output_name, jackknife_theory_covariance=jack_cov, full_theory_covariance=full_cov,
         jackknife_data_covariance=data_cov, shot_noise_rescaling=alpha_best,
         jackknife_theory_precision=jack_prec, full_theory_precision=full_prec,
         N_eff=N_eff_D, full_theory_D_matrix=full_D_est,
         individual_theory_covariances=partial_cov, individual_theory_jackknife_covariances=partial_jack_cov)

print("Saved output covariance matrices as %s"%output_name)
