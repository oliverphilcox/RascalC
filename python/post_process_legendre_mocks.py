## Script to post-process the single-field Legendre binned integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os
from tqdm import trange
from warnings import warn


def post_process_legendre_mocks(mock_cov_file: str, file_root: str, n: int, max_l: int, n_samples: int, outdir: str, skip_r_bins: int = 0, skip_l: int = 0, print_function = print):
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l//2 + 1 # number of multipoles stored

    n_bins0 = n_l * n
    n_bins = (n_l - skip_l) * (n - skip_r_bins)

    mock_cov = np.loadtxt(mock_cov_file) # load external mock covariance matrix
    if mock_cov.shape != (n_bins0, n_bins0): raise ValueError("Unexpected shape of mock covariance matrix")
    mock_cov = mock_cov.reshape(n_l, n, n_l, n) # convert from 2D to 4D with [l, r] ordering in rows and columns
    mock_cov = mock_cov.transpose(1, 0, 3, 2) # change ordering to [r, l] for both rows and columns, like in RascalC results
    mock_cov = mock_cov[skip_r_bins:, :, skip_r_bins:, :] # skip r bins while the dimensions are nicely separated
    if skip_l > 0: mock_cov = mock_cov[:, :-skip_l, :, :-skip_l] # skip multipoles while the dimensions are nicely separated; would be wrong for skip_l=0 so need to skip in this case
    mock_cov = mock_cov.reshape(n_bins, n_bins) # convert back to 2D

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    def load_matrices(index):
        """Load intermediate or full covariance matrices"""
        cov_root = os.path.join(file_root, 'CovMatricesAll/')
        c2 = np.loadtxt(cov_root + 'c2_n%d_l%d_11_%s.txt' % (n, max_l, index))
        c3 = np.loadtxt(cov_root + 'c3_n%d_l%d_1,11_%s.txt' % (n, max_l, index))
        c4 = np.loadtxt(cov_root + 'c4_n%d_l%d_11,11_%s.txt' % (n, max_l, index))

        l_mask = (np.arange(n_l) < n_l - skip_l) # this mask skips last skip_l multipoles
        full_mask = np.append(np.zeros(skip_r_bins * n_l, dtype=bool), np.tile(l_mask, n - skip_r_bins)) # start with zeros and then tile (append to itself n - skip_r_bins times) the l_mask since cov terms are first ordered by r and then by l
        c2, c3, c4 = (a[full_mask][:, full_mask] for a in (c2, c3, c4)) # select rows and columns

        # Now symmetrize and return matrices
        return c2, 0.5*(c3+c3.T), 0.5*(c4+c4.T)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c2f, c3f, c4f = load_matrices('full')

    # Check matrix convergence
    from numpy.linalg import eigvalsh
    eig_c4 = eigvalsh(c4f)
    eig_c2 = eigvalsh(c2f)
    if min(eig_c4) < -1.*min(eig_c2):
        warn("WARNING: 4-point covariance matrix has not converged properly via the eigenvalue test. Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4), min(eig_c2)))

    # Load in partial theoretical matrices
    c2s, c3s, c4s = [], [], []
    for i in trange(n_samples, desc = "Loading full subsamples"):
        c2t, c3t, c4t = load_matrices(i)
        c2s.append(c2t)
        c3s.append(c3t)
        c4s.append(c4t)
    c2s, c3s, c4s = [np.array(a) for a in (c2s, c3s, c4s)]

    # Compute inverted matrix
    def Psi(alpha):
        """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
        c_tot = c2f * alpha**2 + c3f * alpha + c4f
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
        if logdet[0]<0:
            # Remove any dodgy inversions
            return np.inf        
        return np.trace(np.matmul(Psi_alpha, mock_cov))-logdet[1]

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    from scipy.optimize import fmin
    alpha_best = fmin(neg_log_L1, 1.)
    print_function("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)
    alpha = alpha_best # to save editing later

    # Compute full covariance matrices and precision
    full_cov = c4f + c3f*alpha + c2f*alpha**2

    # Check positive definiteness
    if np.any(np.linalg.eigvalsh(full_cov) <= 0): raise ValueError("The full covariance is not positive definite - insufficient convergence")

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    partial_cov = alpha**2 * c2s + alpha * c3s + c4s
    sum_partial_cov = np.sum(partial_cov, axis=0)
    tmp = 0.
    for i in range(n_samples):
        c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
        tmp+=np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
    full_D_est = (n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
    full_prec = np.matmul(np.eye(n_bins)-full_D_est,np.linalg.inv(full_cov))
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    slogdetD = np.linalg.slogdet(full_D_est)
    D_value = slogdetD[0] * np.exp(slogdetD[1] / n_bins)
    if slogdetD[0] < 0:
        print_function("N_eff is negative! Setting to zero")
        N_eff_D = 0.
    else:
        N_eff_D = (n_bins+1.)/D_value+1.
        print_function("Total N_eff Estimate: %.4e" % N_eff_D)

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_Mocks_n%d_l%d.npz' % (n, max_l))
    np.savez(output_name, full_theory_covariance=full_cov,
            shot_noise_rescaling=alpha, full_theory_precision=full_prec,
            N_eff=N_eff_D, full_theory_D_matrix=full_D_est,
            individual_theory_covariances=partial_cov)

    print_function("Saved output covariance matrices as %s" % output_name)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (7, 8, 9):
        print("Usage: python post_process_legendre_mocks.py {MOCK_COV_FILE} {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SKIP_R_BINS} [{SKIP_L}]]")
        sys.exit(1)
            
    mock_cov_file = str(sys.argv[1])
    file_root = str(sys.argv[2])
    n = int(sys.argv[3])
    max_l = int(sys.argv[4])
    n_samples = int(sys.argv[5])
    outdir = str(sys.argv[6])
    from utils import get_arg_safe
    skip_r_bins = get_arg_safe(7, int, 0)
    skip_l = get_arg_safe(8, int, 0)

    post_process_legendre_mocks(mock_cov_file, file_root, n, max_l, n_samples, outdir, skip_r_bins, skip_l)