## Script to post-process the single-field Legendre binned integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys, os
from tqdm import trange


def post_process_3pcf(file_root: str, n: int, max_l: int, n_samples: int, outdir: str, alpha: float = 1, print_function = print):
    m = max_l+1
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    def symmetrize(mat):
        """ Add in symmetries to matrices """
        out_mat = np.zeros_like(mat)
        for i in range(len(mat)//m):
            a = i//n
            b = i%n
            for j in range(len(mat)//m):
                c = j//n
                d = j%n
                # Add to all relevant bins
                these_mat = mat[i*m:(i+1)*m,j*m:(j+1)*m]*0.25
                out_mat[(a*n+b)*m:(a*n+b)*m+m,(c*n+d)*m:(c*n+d)*m+m]+=these_mat
                out_mat[(b*n+a)*m:(b*n+a)*m+m,(c*n+d)*m:(c*n+d)*m+m]+=these_mat
                out_mat[(b*n+a)*m:(b*n+a)*m+m,(d*n+c)*m:(d*n+c)*m+m]+=these_mat
                out_mat[(a*n+b)*m:(a*n+b)*m+m,(d*n+c)*m:(d*n+c)*m+m]+=these_mat
        return 0.5*(out_mat+out_mat.T)

    def load_matrices(index):
        """Load intermediate or full covariance matrices"""
        cov_root = os.path.join(file_root, '3PCFCovMatricesAll/')
        c3_0 = np.loadtxt(cov_root+'c3_n%d_l%d_0_%s.txt'%(n,max_l,index))
        c4_0 = np.loadtxt(cov_root+'c4_n%d_l%d_0_%s.txt'%(n,max_l,index))
        c5_0 = np.loadtxt(cov_root+'c5_n%d_l%d_0_%s.txt'%(n,max_l,index))
        c6_0 = np.loadtxt(cov_root+'c6_n%d_l%d_0_%s.txt'%(n,max_l,index))
        c3_1 = np.loadtxt(cov_root+'c3_n%d_l%d_1_%s.txt'%(n,max_l,index))
        c4_1 = np.loadtxt(cov_root+'c4_n%d_l%d_1_%s.txt'%(n,max_l,index))
        c5_1 = np.loadtxt(cov_root+'c5_n%d_l%d_1_%s.txt'%(n,max_l,index))
        c6_1 = np.loadtxt(cov_root+'c6_n%d_l%d_1_%s.txt'%(n,max_l,index))
        
        c3 = c3_0+c3_1
        c4 = c4_0+c4_1
        c5 = c5_0+c5_1
        c6 = c6_1 # do not include c6_0 term here
        
        # Now symmetrize and return matrices
        return symmetrize(c3),symmetrize(c4),symmetrize(c5),symmetrize(c6)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c3,c4,c5,c6=load_matrices('full')

    # Compute full covariance matrices and precision
    full_cov = c6+c5*alpha+c4*alpha**2.+c3*alpha**3.
    n_bins = len(c6)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial theoretical matrices
    c3s, c4s, c5s, c6s = [], [], [], []
    for i in trange(n_samples, desc="Loading full subsamples"):
        c3t, c4t, c5t, c6t = load_matrices(i)
        c3s.append(c3t)
        c4s.append(c4t)
        c5s.append(c5t)
        c6s.append(c6t)
    c3s, c4s, c5s, c6s = [np.array(a) for a in (c3s, c4s, c5s, c6s)]
        
    partial_cov = alpha**3 * c3s + alpha**2 * c4s + alpha * c5s + c6s
    sum_partial_cov = np.sum(partial_cov, axis=0)

    tmp=0.

    for _ in range(1):
        for i in range(n_samples):
            c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
            try:
                tmp+=np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
            except np.linalg.linalg.LinAlgError:
                print_function("Could not invert submatrix, so setting overall precision to zero. Matrix is not fully converged")
                tmp = np.inf
        if tmp==np.inf:
            # Couldn't estimate precision so break loop
            full_prec = np.zeros_like(full_cov)
            full_D_est = np.zeros_like(full_cov)
            N_eff_D = 0.
            continue
                
        full_D_est=(n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
        full_prec = np.matmul(np.eye(n_bins)-full_D_est,np.linalg.inv(full_cov))
        print_function("Full precision matrix estimate computed")    

        # Now compute effective N:
        slogdetD=np.linalg.slogdet(full_D_est)
        D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
        if slogdetD[0]<0:
            print_function("N_eff is negative! Setting to zero")
            N_eff_D = 0.
        else:
            N_eff_D = (n_bins+1.)/D_value+1.
            print_function("Total N_eff Estimate: %.4e"%N_eff_D)        

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_3PCF_n%d_l%d.npz'%(n,max_l))
    np.savez(output_name,full_theory_covariance=full_cov,
            shot_noise_rescaling=alpha,full_theory_precision=full_prec,
            N_eff=N_eff_D,full_theory_D_matrix=full_D_est,
            individual_theory_covariances=partial_cov)

    print_function("Saved output covariance matrices as %s"%output_name)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (6, 7):
        print("Usage: python post_process_3pcf.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]")
        sys.exit(1)
            
    file_root = str(sys.argv[1])
    n = int(sys.argv[2])
    max_l = int(sys.argv[3])
    n_samples = int(sys.argv[4])
    outdir = str(sys.argv[5])
    from .utils import get_arg_safe
    alpha = get_arg_safe(6, float, 1)

    post_process_3pcf(file_root, n, max_l, n_samples, outdir, alpha)