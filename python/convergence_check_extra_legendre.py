## Script to perform an extra convergence check on full integrals in Legendre mode
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Determines single-field vs multi-field automatically

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (5, 6, 7, 8): # if too few or parity is wrong
    print("Usage: python convergence_check_extra.py {N_R_BINS} {MAX_L} {COVARIANCE_INPUT_DIR} {N_SUBSAMPLES} [{ALPHA} [{SKIP_R_BINS} [{SKIP_L}]]]")
    sys.exit()

n = int(sys.argv[1])
max_l = int(sys.argv[2])
input_root = str(sys.argv[3])
n_samples = int(sys.argv[4])
alpha = float(sys.argv[5]) if len(sys.argv) >= 6 else 1
skip_bins = int(sys.argv[6]) if len(sys.argv) >= 7 else 0
skip_l = int(sys.argv[7]) if len(sys.argv) >= 8 else 0

# mask for skipping bins
r_mask = np.arange(n) >= skip_bins

input_root_all = os.path.join(input_root, 'CovMatricesAll/')

# methods to assess similarity
def rms_eig_inv_test_covs(C1, C2):
    Psi1 = np.linalg.inv(C1)
    N = len(C2)
    tmp = Psi1.dot(C2) - np.eye(N)
    return np.sqrt(np.sum(np.diag(tmp.dot(tmp))) / N)

def KL_div_covs(C1, C2):
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return (np.trace(Psi1C2) - len(C2) - np.log(np.linalg.det(Psi1C2)))/2

def symmetrized(A): # symmetrize a 2D array
    return 0.5 * (A + A.T)

# input indices
I1 = [1,1,1,1,1,2,2]
I2 = [1,2,2,2,1,1,2]
I3 = [1,1,2,1,2,2,2]
I4 = [1,1,1,2,2,2,2]

for ii in range(len(I1)): # loop over all field combinations
    index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
    index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
    index2="%d%d"%(I1[ii],I2[ii])

    # full integrals
    c2, c3, c4 = [], [], []
    # read
    for i in range(n_samples):
        try:
            c2.append(np.diag(np.loadtxt(input_root_all+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i))))
        except (FileNotFoundError, IOError): break # end loop if c2 full not found
        c3.append(np.loadtxt(input_root_all+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i)))
        c4.append(np.loadtxt(input_root_all+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i)))
    if len(c2) == 0: break # end loop if no full integral has been found
    if len(c2) < n_samples:
        print("Some %s full samples missing: expected %d, found %d" % (index4, n_samples_tot, len(c2)))
        break # end loop like above
    c2, c3, c4 = np.array(c2), np.array(c3), np.array(c4)
    # construct boolean mask
    N = np.shape(c2)[1]
    assert N % n == 0, "Number of bins mismatch"
    nrepeat = N // n - skip_l
    full_mask = np.append(np.repeat(r_mask, nrepeat), np.zeros(n * skip_l, dtype=bool)) # repeat the r_mask and append zeros since cov terms are first ordered by l
    # construct averages in halves
    c2_first, c3_first, c4_first = [np.mean(a[:n_samples//2, full_mask][:, :, full_mask], axis=0) for a in (c2, c3, c4)]
    cov_first = c2_first * alpha**2 + symmetrized(c3_first) * alpha + symmetrized(c4_first)
    c2_second, c3_second, c4_second = [np.mean(a[n_samples//2:, full_mask][:, :, full_mask], axis=0) for a in (c2, c3, c4)]
    cov_second = c2_second * alpha**2 + symmetrized(c3_second) * alpha + symmetrized(c4_second)
    print("%s full: RMS eigenvalues of inverse tests for cov half-estimates are %.2e and %.2e" % (index4, rms_eig_inv_test_covs(cov_first, cov_second), rms_eig_inv_test_covs(cov_second, cov_first)))
    print("%s full: KL divergences between cov half-estimates are %.2e and %.2e" % (index4, KL_div_covs(cov_first, cov_second), KL_div_covs(cov_second, cov_first)))
