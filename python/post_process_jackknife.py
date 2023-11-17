## Script to post-process the single-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import numpy as np
import sys,os
from tqdm import trange
from warnings import warn


def post_process_jackknife(jackknife_file: str, weight_dir: str, file_root: str, m: int, n_samples: int, outdir: str, skip_r_bins: int = 0, print_function = print):
    skip_bins = skip_r_bins * m

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load jackknife xi estimates from data
    print_function("Loading correlation function jackknife estimates from %s"%jackknife_file)
    xi_jack = np.loadtxt(jackknife_file, skiprows=2)[:, skip_bins:]
    n_bins = xi_jack.shape[1] # total bins
    n_jack = xi_jack.shape[0] # total jackknives
    n = n_bins // m + skip_r_bins # radial bins

    weight_file = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack))
    RR_file = os.path.join(weight_dir, 'binned_pair_counts_n%d_m%d_j%d_11.dat'%(n,m,n_jack))

    print_function("Loading weights file from %s"%weight_file)
    weights = np.loadtxt(weight_file)[:, 1+skip_bins:]

    # First exclude any dodgy jackknife regions
    good_jk = np.where(np.all(np.isfinite(xi_jack), axis=1))[0] # all xi in jackknife have to be normal numbers
    print_function("Using %d out of %d jackknives"%(len(good_jk),n_jack))
    xi_jack = xi_jack[good_jk]
    weights = weights[good_jk]
    weights /= np.sum(weights, axis=0) # renormalize weights after possibly discarding some jackknives

    # Compute data covariance matrix
    print_function("Computing data covariance matrix")
    mean_xi = np.sum(xi_jack*weights,axis=0)
    tmp = weights*(xi_jack-mean_xi)
    data_cov = np.matmul(tmp.T,tmp)
    denom = np.matmul(weights.T,weights)
    data_cov /= (np.ones_like(denom)-denom)

    print_function("Loading weights file from %s"%RR_file)
    RR=np.loadtxt(RR_file)

    def load_matrices(index,jack=True):
        """Load intermediate or full covariance matrices"""
        if jack:
            cov_root = os.path.join(file_root, 'CovMatricesJack/')
        else:
            cov_root = os.path.join(file_root, 'CovMatricesAll/')
        c2 = np.diag(np.loadtxt(cov_root+'c2_n%d_m%d_11_%s.txt'%(n,m,index))[skip_bins:])
        c3 = np.loadtxt(cov_root+'c3_n%d_m%d_1,11_%s.txt'%(n,m,index))[skip_bins:, skip_bins:]
        c4 = np.loadtxt(cov_root+'c4_n%d_m%d_11,11_%s.txt'%(n,m,index))[skip_bins:, skip_bins:]
        if jack:
            EEaA1 = np.loadtxt(cov_root+'EE1_n%d_m%d_11_%s.txt' %(n,m,index))[:, skip_bins:]
            EEaA2 = np.loadtxt(cov_root+'EE2_n%d_m%d_11_%s.txt' %(n,m,index))[:, skip_bins:]
            RRaA1 = np.loadtxt(cov_root+'RR1_n%d_m%d_11_%s.txt' %(n,m,index))[:, skip_bins:]
            RRaA2 = np.loadtxt(cov_root+'RR2_n%d_m%d_11_%s.txt' %(n,m,index))[:, skip_bins:]
        
            # Compute disconnected term
            w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
            w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
            diff1 = EEaA1-w_aA1*EEaA1.sum(axis=0)
            diff2 = EEaA2-w_aA2*EEaA2.sum(axis=0)
            RRaRRb = np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
            fact = np.ones_like(c4)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
            cx = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
            c4+=cx
        
        # Now symmetrize and return matrices
        return c2,0.5*(c3+c3.T),0.5*(c4+c4.T)

    # Load in full jackknife theoretical matrices
    print_function("Loading best estimate of jackknife covariance matrix")
    c2j,c3j,c4j=load_matrices('full')

    # Check matrix convergence
    from numpy.linalg import eigvalsh
    eig_c4 = eigvalsh(c4j)
    eig_c2 = eigvalsh(c2j)
    if min(eig_c4)<-1.*min(eig_c2):
        warn("Jackknife 4-point covariance matrix has not converged properly via the eigenvalue test. Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4), min(eig_c2)))

    # Load in partial jackknife theoretical matrices
    c2s, c3s, c4s = [], [], []
    for i in trange(n_samples, desc="Loading jackknife subsamples"):
        c2, c3, c4 = load_matrices(i)
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
        tmp=0.
        for i in range(n_samples):
            c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
            tmp+=np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
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
    print_function("Optimizing for the shot-noise rescaling parameter")
    from scipy.optimize import fmin
    alpha_best = fmin(neg_log_L1,1.)
    print_function("Optimization complete - optimal rescaling parameter is %.6f"%alpha_best)

    # Compute jackknife and full covariance matrices
    jack_cov = c4j+c3j*alpha_best+c2j*alpha_best**2.
    jack_prec = Psi(alpha_best)
    c2f,c3f,c4f=load_matrices('full',jack=False)
    full_cov = c4f+c3f*alpha_best+c2f*alpha_best**2.

    # Check positive definiteness
    if np.any(np.linalg.eigvalsh(full_cov) <= 0): raise ValueError("The full covariance is not positive definite - insufficient convergence")

    # Check convergence
    eig_c4f = eigvalsh(c4f)
    eig_c2f = eigvalsh(c2f)
    if min(eig_c4f)<min(eig_c2f)*-1.:
        raise ValueError("Full 4-point covariance matrix has not converged properly via the eigenvalue test. Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4f), min(eig_c2f)))

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial jackknife theoretical matrices
    c2fs, c3fs, c4fs = [], [], []
    for i in trange(n_samples, desc="Loading full subsamples"):
        c2, c3, c4 = load_matrices(i, jack=False)
        c2fs.append(c2)
        c3fs.append(c3)
        c4fs.append(c4)
    c2fs, c3fs, c4fs = [np.array(a) for a in (c2fs, c3fs, c4fs)]
    partial_cov = alpha_best**2 * c2fs + alpha_best * c3fs + c4fs
    sum_partial_cov = np.sum(partial_cov, axis=0)
    tmp=0.
    for i in range(n_samples):
        c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
        tmp+=np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
    full_D_est=(n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
    full_prec = np.matmul(np.eye(n_bins)-full_D_est,np.linalg.inv(full_cov))
    print_function("Full precision matrix estimate computed")    

    # Now compute effective N:
    slogdetD = np.linalg.slogdet(full_D_est)
    D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
    N_eff_D = (n_bins+1.)/D_value+1.
    print_function("Total N_eff Estimate: %.4e" % N_eff_D)

    # Jackknife covariance for posterity
    partial_jack_cov = alpha_best**2 * c2s + alpha_best * c3s + c4s

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Jackknife_n%d_m%d_j%d.npz' % (n, m, n_jack))
    np.savez(output_name, jackknife_theory_covariance=jack_cov, full_theory_covariance=full_cov,
            jackknife_data_covariance=data_cov, shot_noise_rescaling=alpha_best,
            jackknife_theory_precision=jack_prec, full_theory_precision=full_prec,
            N_eff=N_eff_D, full_theory_D_matrix=full_D_est,
            individual_theory_covariances=partial_cov, individual_theory_jackknife_covariances=partial_jack_cov)

    print_function("Saved output covariance matrices as %s"%output_name)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (7, 8):
        print("Usage: python post_process_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SKIP_R_BINS}]")
        sys.exit(1)

    jackknife_file = str(sys.argv[1])
    weight_dir = str(sys.argv[2])
    file_root = str(sys.argv[3])
    m = int(sys.argv[4])
    n_samples = int(sys.argv[5])
    outdir = str(sys.argv[6])
    from utils import get_arg_safe
    skip_r_bins = get_arg_safe(7, int, 0)

    post_process_jackknife(jackknife_file, weight_dir, file_root, m, n_samples, outdir, skip_r_bins)