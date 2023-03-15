## Script to post-process the multi-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (6, 8):
    print("Usage: python post_process_default_multi.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2}]")
    sys.exit(1)

file_root = str(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
n_samples = int(sys.argv[4])
outdir = str(sys.argv[5])
alpha_1 = float(sys.argv[6]) if len(sys.argv) >= 7 else 1
alpha_2 = float(sys.argv[7]) if len(sys.argv) >= 8 else 1
n_bins = n * m

alphas = [alpha_1, alpha_2]

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

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
                print("Reading in integral components for C_{%s}, which used %.2e pairs, %.2e triples and %.2e quads of particles" % (index4,total_counts[0], total_counts[1], total_counts[2]))
            except (FileNotFoundError, IOError): pass
        else:
            pass
            #print("Reading in integral components for C_{%s}, iteration %s"%(index4,suffix))

        # Load full integrals
        c2 = np.diag(np.loadtxt(file_root_all + 'c2_n%d_m%d_%s_%s.txt' % (n, m, index2, suffix)))
        c3 = np.loadtxt(file_root_all + 'c3_n%d_m%d_%s_%s.txt' % (n, m, index3, suffix))
        c4 = np.loadtxt(file_root_all + 'c4_n%d_m%d_%s_%s.txt' % (n, m, index4, suffix))

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

    c_tot = np.zeros([3,3,n*m,n*m]) # array with each individual covariance accessible
    c_comb = np.zeros([3*n*m,3*n*m]) # full array suitable for inversion

    for j1 in range(3):
        ind1,ind2 = cov_indices[j1]
        alpha1, alpha2 = alphas[ind1], alphas[ind2]
        for j2 in range(3):
            ind3,ind4 = cov_indices[j2]
            tmp = construct_fields(ind1, ind2, ind3, ind4, alpha1, alpha2)
            c_tot[j1,j2] = tmp
            c_comb[j1*n*m:(j1+1)*n*m,j2*n*m:(j2+1)*n*m] = tmp

    return c_tot,0.5*(c_comb+c_comb.T) # add all remaining symmetries

# Load full matrices
c_tot, c_comb = matrix_readin()
n_bins = len(c_tot[0,0])

# Load subsampled matrices (all submatrices combined)
c_subsamples=[]
for i in range(n_samples):
    _,tmp=matrix_readin(i)
    c_subsamples.append(tmp)
c_subsamples = np.array(c_subsamples)

# Now compute all precision matrices
iden = np.eye(len(c_comb))

N_eff = np.zeros([2,2,2,2])
D_est = np.zeros_like(c_tot)

def compute_precision(entire_matrix, subsamples):
    summ=0.
    sum_subsamples = np.sum(subsamples, axis=0)
    for i in range(n_samples):
        c_excl_i = (sum_subsamples - subsamples[i]) / (n_samples - 1)
        summ+=np.matmul(np.linalg.inv(c_excl_i), subsamples[i])
    D_est = (summ/n_samples-iden)*(n_samples-1.)/n_samples
    logdetD = np.linalg.slogdet(D_est)
    if logdetD[0]<0:
        N_eff_D = 0.
    else:
        D_value = logdetD[0]*np.exp(logdetD[1]/n_bins)
        N_eff_D = (n_bins+1.)/D_value+1.
    precision = np.matmul(iden-D_est,np.linalg.inv(entire_matrix))
    return precision,N_eff_D,D_est

print("Computing precision matrices and N_eff")
prec_comb,N_eff,D_est = compute_precision(c_comb, c_subsamples)

output_name =os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Default_n%d_m%d.npz'%(n,m))
np.savez(output_name,full_theory_covariance=c_comb,
         all_covariances = c_tot,
         shot_noise_rescaling=[alpha_1,alpha_2],
         full_theory_precision=prec_comb,
         N_eff=N_eff, full_theory_D_matrix = D_est,
         individual_theory_covariances = c_subsamples)

print("Saved output covariance matrices as %s"%output_name)
