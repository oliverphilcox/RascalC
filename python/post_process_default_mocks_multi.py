## Script to post-process the multi-field integrals computed by the C++ code. This computes two shot-noise rescaling parameters, alphas, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os
from tqdm import trange


def post_process_default_mocks_multi(mock_cov_file: str, file_root: str, n: int, m: int, n_samples: int, outdir: str, skip_r_bins: int = 0, print_function = print):
    skip_bins = skip_r_bins * m
    n_bins = n * m - skip_bins

    skip_mask = np.tile(np.arange(n * m) >= skip_bins, 3) # the mask gives False for first skip_bins, all this repeating 3 times; 3 is number of correlations for 2 tracers
    mock_cov = np.loadtxt(mock_cov_file)[skip_mask][:, skip_mask] # load external mock covariance matrix, select bins on both axes

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    def load_matrices(index, field):
        """Load intermediate or full autocovariance matrices.
        The field parameter controls which field covariance matrix to load"""
        cov_root = os.path.join(file_root, 'CovMatricesAll/')
        suffix2 = '_n%d_m%d_%s%s_%s.txt' % (n, m, field, field, index)
        suffix3 = '_n%d_m%d_%s,%s%s_%s.txt' % (n, m, field, field, field, index)
        suffix4 = '_n%d_m%d_%s%s,%s%s_%s.txt' % (n, m, field, field, field, field, index)
        c2 = np.diag(np.loadtxt(cov_root + 'c2' + suffix2)[skip_bins:])
        c3 = np.loadtxt(cov_root + 'c3' + suffix3)[skip_bins:, skip_bins:]
        c4 = np.loadtxt(cov_root + 'c4' + suffix4)[skip_bins:, skip_bins:]

        # Now symmetrize and return matrices
        return c2, 0.5*(c3+c3.T), 0.5*(c4+c4.T)

    # Load autocovariance from mock covariance
    mock_cov_11 = mock_cov[:n_bins, :n_bins]
    mock_cov_22 = mock_cov[2*n_bins:, 2*n_bins:]
    auto_mock_cov = [mock_cov_11, mock_cov_22]

    alpha_best = np.zeros(2)

    indices = ['11', '22']

    ## Optimize for alpha_1 and alpha_2 separately.
    for i, index in enumerate(indices):

        # Load mock auto-covariance matrix
        this_mock_cov = auto_mock_cov[i]

        # Load in full jackknife theoretical matrices
        print_function("Loading best estimate of covariance matrix for field %d" % (i+1))
        c2f, c3f, c4f = load_matrices('full', i+1)

        # Check matrix convergence
        from numpy.linalg import eigvalsh
        eig_c4 = eigvalsh(c4f)
        eig_c2 = eigvalsh(c2f)
        if min(eig_c4) < -1.*min(eig_c2):
            raise ValueError("4-point covariance matrix has not converged properly via the eigenvalue test. Min eigenvalue of C4 = %.2e, min eigenvalue of C2 = %.2e" % (min(eig_c4), min(eig_c2)))

        # Load in partial jackknife theoretical matrices
        c2s, c3s, c4s = [], [], []
        for j in trange(n_samples, desc=f"Loading field {i+1} subsamples"):
            c2, c3, c4 = load_matrices(j, i+1)
            c2s.append(c2)
            c3s.append(c3)
            c4s.append(c4)
        c2s, c3s, c4s = [np.array(a) for a in (c2s, c3s, c4s)]

        # Compute inverted matrix
        def Psi(alpha):
            """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
            c_tot = c2f * alpha**2. + c3f * alpha + c4f
            partial_cov = alpha**2 * c2s + alpha * c3s + c4s
            sum_partial_cov = np.sum(partial_cov, axis=0)
            tmp = 0.
            for i in range(n_samples):
                c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
                tmp += np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
            D_est = (n_samples - 1) / n_samples * (-np.eye(n_bins) + tmp / n_samples)
            Psi = np.matmul(np.eye(n_bins) - D_est, np.linalg.inv(c_tot))
            return Psi

        def neg_log_L1(alpha):
            """Return negative log L1 likelihood between data and theory covariance matrices"""
            Psi_alpha = Psi(alpha)
            logdet = np.linalg.slogdet(Psi_alpha)
            if logdet[0] < 0:
                # Remove any dodgy inversions
                return np.inf
            return np.trace(np.matmul(Psi_alpha, this_mock_cov)) - logdet[1]

        # Now optimize for shot-noise rescaling parameter alpha
        print_function("Optimizing for the shot-noise rescaling parameter alpha_%d" % (i+1))
        from scipy.optimize import fmin
        optimal_alpha = fmin(neg_log_L1, 1.)
        print_function("Optimization complete for field %d - optimal rescaling parameter is alpha_%d = %.6f" % (i+1, i+1, optimal_alpha))

        alpha_best[i] = optimal_alpha

    # input indices
    I1 = [1,1,1,1,1,2,2]
    I2 = [1,2,2,2,1,1,2]
    I3 = [1,1,2,1,2,2,2]
    I4 = [1,1,1,2,2,2,2]

    def matrix_readin(suffix='full'):
        """Read in multi-field covariance matrices. This returns lists of covariance matrices and a combined covariance matrix."""

        ## Define arrays for covariance matrices
        c2s = np.zeros([2, 2, n_bins, n_bins])
        c3s = np.zeros([2, 2, 2, n_bins, n_bins])
        c4s = np.zeros([2, 2, 2, 2, n_bins, n_bins])
        ## Normalization arrays for covariance matrices
        n2s = np.zeros([2, 2])
        n3s = np.zeros([2, 2, 2])
        n4s = np.zeros([2, 2, 2, 2])

        for ii in range(len(I1)):
            index4 = "%d%d,%d%d" % (I1[ii], I2[ii], I3[ii], I4[ii])
            index3 = "%d,%d%d" % (I2[ii], I1[ii], I3[ii])
            index2 = "%d%d" % (I1[ii], I2[ii])

            j1, j2, j3,j4=I1[ii]-1, I2[ii]-1, I3[ii]-1, I4[ii]-1 # internal indexing

            # Define input files
            file_root_all = os.path.join(file_root, 'CovMatricesAll/')

            if suffix=='full':
                counts_file = file_root_all + 'total_counts_n%d_m%d_%s.txt' % (n, m, index4)
                # Load total number of counts
                try:
                    total_counts=np.loadtxt(counts_file)
                    print_function("Reading in integral components for C_{%s}, which used %.2e pairs, %.2e triples and %.2e quads of particles" % (index4,total_counts[0], total_counts[1], total_counts[2]))
                except (FileNotFoundError, IOError): pass
            else:
                pass
                #print_function("Reading in integral components for C_{%s}, iteration %s"%(index4,suffix))

            # Load full integrals
            c2 = np.diag(np.loadtxt(file_root_all + 'c2_n%d_m%d_%s_%s.txt' % (n, m, index2, suffix))[skip_bins:])
            c3 = np.loadtxt(file_root_all + 'c3_n%d_m%d_%s_%s.txt' % (n, m, index3, suffix))[skip_bins:, skip_bins:]
            c4 = np.loadtxt(file_root_all + 'c4_n%d_m%d_%s_%s.txt' % (n, m, index4, suffix))[skip_bins:, skip_bins:]

            # Now save components
            c2s[j1, j2] += c2
            n2s[j1, j2] += 1
            c3s[j2, j1, j3] += c3
            n3s[j2, j1, j3] += 1
            # will deal with c4s/n4s later

            # c2 symmetry - indices interchanged, ensures matrix symmetry if they are equal
            c2s[j2, j1] += c2
            n2s[j2, j1] += 1

            # c3 symmetry - last two indices interchanged, ensures matrix symmetry if they are equal
            c3s[j2, j3, j1] += c3.T
            n3s[j2, j3, j1] += 1
            
            # All symmetries possible for c4 without transpositions
            permutations4 = ((j1, j2, j3, j4), # original
                            (j2, j1, j3, j4), # first two indices interchanged
                            (j1, j2, j4, j3), # last two indices interchanged
                            (j2, j1, j4, j3), # first and last two indices interchanged at the same time
                            )
            
            for (i1, i2, i3, i4) in permutations4:
                c4s[i1, i2, i3, i4] += c4
                n4s[i1, i2, i3, i4] += 1
                # now swap indices and transpose
                c4s[i3, i4, i1, i2] += c4.T
                n4s[i3, i4, i1, i2] += 1
        
        # normalize the covariances
        c2s /= n2s[:, :, None, None]
        c3s /= n3s[:, :, :, None, None]
        c4s /= n4s[:, :, :, :, None, None]

        def construct_fields(j1, j2, j3, j4, alpha1, alpha2):
            # Reconstruct the full field for given input fields and rescaling parameters

            # Create kronecker deltas
            d_xw = (j1 == j4)
            d_xz = (j1 == j3)
            d_yw = (j2 == j4)
            d_yz = (j2 == j3)

            full = c4s[j1, j2, j3, j4] + 0.25 * alpha1 * (d_xw * c3s[j1, j2, j3] + d_xz * c3s[j1, j2, j4]) + 0.25 * alpha2 * (d_yw * c3s[j2, j1, j3] + d_yz * c3s[j2, j1, j4]) + 0.5 * alpha1 * alpha2 * (d_xw * d_yz + d_xz * d_yw) * c2s[j1, j2]
            return full

        # Index in ordering (P_11,P_12,P_22)
        cov_indices = [[0, 0], [0, 1], [1, 1]]

        c_tot = np.zeros([3, 3, n_bins, n_bins]) # array with each individual covariance accessible
        c_comb = np.zeros([3*n_bins, 3*n_bins]) # full array suitable for inversion

        for j1 in range(3):
            ind1, ind2 = cov_indices[j1]
            alpha1, alpha2 = alpha_best[[ind1, ind2]]
            for j2 in range(3):
                ind3,ind4 = cov_indices[j2]
                tmp = construct_fields(ind1, ind2, ind3, ind4, alpha1, alpha2)
                c_tot[j1, j2] = tmp
                c_comb[j1*n_bins:(j1+1)*n_bins, j2*n_bins:(j2+1)*n_bins] = tmp

        return c_tot, 0.5*(c_comb+c_comb.T) # add all remaining symmetries

    # Load full matrices
    c_tot, c_comb = matrix_readin()
    n_bins = len(c_tot[0,0])

    # Load subsampled matrices (all submatrices combined)
    c_subsamples=[]
    for i in trange(n_samples, desc="Loading full subsamples"):
        _, tmp = matrix_readin(i)
        c_subsamples.append(tmp)
    c_subsamples = np.array(c_subsamples)

    # Now compute all precision matrices
    iden = np.eye(len(c_comb))

    N_eff = np.zeros([2, 2, 2, 2])
    D_est = np.zeros_like(c_tot)

    def compute_precision(entire_matrix, subsamples):
        summ=0.
        sum_subsamples = np.sum(subsamples, axis=0)
        for i in range(n_samples):
            c_excl_i = (sum_subsamples - subsamples[i]) / (n_samples - 1)
            summ+=np.matmul(np.linalg.inv(c_excl_i), subsamples[i])
        D_est = (summ / n_samples - iden) * (n_samples-1.) / n_samples
        logdetD = np.linalg.slogdet(D_est)
        if logdetD[0]<0:
            N_eff_D = 0.
        else:
            D_value = logdetD[0] * np.exp(logdetD[1] / n_bins)
            N_eff_D = (n_bins+1.) / D_value + 1.
        precision = np.matmul(iden-D_est, np.linalg.inv(entire_matrix))
        return precision,N_eff_D,D_est

    print_function("Computing precision matrices and N_eff")
    prec_comb,N_eff,D_est = compute_precision(c_comb, c_subsamples)

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Default_Mocks_n%d_m%d.npz' % (n, m))
    np.savez(output_name, full_theory_covariance=c_comb,
            all_covariances = c_tot,
            shot_noise_rescaling=alpha_best,
            full_theory_precision=prec_comb,
            N_eff=N_eff, full_theory_D_matrix = D_est,
            individual_theory_covariances = c_subsamples)

    print_function("Saved output covariance matrices as %s" % output_name)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (7, 8):
        print("Usage: python post_process_default_mocks_multi.py {MOCK_COV_FILE} {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SKIP_R_BINS}]")
        sys.exit(1)
            
    mock_cov_file = str(sys.argv[1])
    file_root = str(sys.argv[2])
    n = int(sys.argv[3])
    m = int(sys.argv[4])
    n_samples = int(sys.argv[5])
    outdir = str(sys.argv[6])
    from utils import get_arg_safe
    skip_r_bins = get_arg_safe(7, int, 0) # convert from radial to total number of bins right away

    post_process_default_mocks_multi(mock_cov_file, file_root, n, m, n_samples, outdir, skip_r_bins)