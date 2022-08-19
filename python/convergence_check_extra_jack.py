## Script to perform an extra convergence check on jackknife integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Determines single-field vs multi-field automatically

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (7, 8, 9): # if too few or parity is wrong
    print("Usage: python convergence_check_extra_jack.py {N_R_BINS} {N_MU_BINS} {COVARIANCE_INPUT_DIR} {N_SUBSAMPLES} {JACKKNIFE_WEIGHTS_FILE} {RR_FILE} [{ALPHA} [{SKIP_R_BINS}]]")
    sys.exit()

n = int(sys.argv[1])
m = int(sys.argv[2])
input_root = str(sys.argv[3])
n_samples = int(sys.argv[4])
weight_file = str(sys.argv[5])
RR_file = str(sys.argv[6])
alpha = float(sys.argv[7]) if len(sys.argv) >= 8 else 1
skip_bins = int(sys.argv[8]) * m if len(sys.argv) >= 9 else 0

input_root_jack = os.path.join(input_root, 'CovMatricesJack/')

weights = np.loadtxt(weight_file)[:,1:]
RR = np.loadtxt(RR_file)

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

# input indices
I1 = [1,1,1,1,1,2,2]
I2 = [1,2,2,2,1,1,2]
I3 = [1,1,2,1,2,2,2]
I4 = [1,1,1,2,2,2,2]

for ii in range(len(I1)): # loop over all field combinations
    index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
    index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
    index2="%d%d"%(I1[ii],I2[ii])

    # jackknife integrals
    c2j, c3j, c4j = [], [], []
    EEaA1, EEaA2, RRaA1, RRaA2 = [], [], [], []
    # read
    for i in range(n_samples):
        try:
            c2j.append(np.diag(np.loadtxt(input_root_jack+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i))))
        except (FileNotFoundError, IOError): break # end loop if c2 jack not found
        c3j.append(np.loadtxt(input_root_jack+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i)))
        c4j.append(np.loadtxt(input_root_jack+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i)))
        # cxj components
        EEaA1.append(np.loadtxt(input_root_jack+'EE1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        EEaA2.append(np.loadtxt(input_root_jack+'EE2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        RRaA1.append(np.loadtxt(input_root_jack+'RR1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        RRaA2.append(np.loadtxt(input_root_jack+'RR2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
    if len(c2j) == 0: break # terminate the loop if no jack integral has been found
    if len(c2j) < n_samples:
        print("Some %s jack samples missing: expected %d, found %d" % (index4, n_samples_tot, len(c2j)))
        continue # skip the rest of the loop like above
    c2j, c3j, c4j = np.array(c2j), np.array(c3j), np.array(c4j)
    EEaA1, EEaA2, RRaA1, RRaA2 = np.array(EEaA1), np.array(EEaA2), np.array(RRaA1), np.array(RRaA2)
    # construct averages in halves
    c2_first, c3_first, c4_first = [np.mean(a[:n_samples//2, skip_bins:, skip_bins:], axis=0) for a in (c2j, c3j, c4j)]
    EEaA1_first, EEaA2_first, RRaA1_first, RRaA2_first = [np.mean(a[:n_samples//2, skip_bins:], axis=0) for a in (EEaA1, EEaA2, RRaA1, RRaA2)]
    w_aA1 = RRaA1_first/np.sum(RRaA1_first,axis=0)
    w_aA2 = RRaA2_first/np.sum(RRaA2_first,axis=0)
    diff1 = EEaA1_first-w_aA1*EEaA1_first.sum(axis=0)
    diff2 = EEaA2_first-w_aA2*EEaA2_first.sum(axis=0)
    RRaRRb = np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
    fact = np.ones_like(c4_first)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
    cx = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
    c4_first+=cx
    cov_first = c2_first * alpha**2 + c3_first * alpha + c4_first
    # second half
    c2_second, c3_second, c4_second = [np.mean(a[n_samples//2:, skip_bins:, skip_bins:], axis=0) for a in (c2j, c3j, c4j)]
    EEaA1_second, EEaA2_second, RRaA1_second, RRaA2_second = [np.mean(a[n_samples//2:, skip_bins:], axis=0) for a in (EEaA1, EEaA2, RRaA1, RRaA2)]
    w_aA1 = RRaA1_second/np.sum(RRaA1_second,axis=0)
    w_aA2 = RRaA2_second/np.sum(RRaA2_second,axis=0)
    diff1 = EEaA1_second-w_aA1*EEaA1_second.sum(axis=0)
    diff2 = EEaA2_second-w_aA2*EEaA2_second.sum(axis=0)
    RRaRRb = np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
    fact = np.ones_like(c4_second)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
    cx = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
    c4_second+=cx
    cov_second = c2_second * alpha**2 + c3_second * alpha + c4_second
    print("%s jack: RMS eigenvalues of inverse tests for cov half-estimates are %.2e and %.2e" % (index4, rms_eig_inv_test_covs(cov_first, cov_second), rms_eig_inv_test_covs(cov_second, cov_first)))
    print("%s jack: KL divergences between cov half-estimates are %.2e and %.2e" % (index4, KL_div_covs(cov_first, cov_second), KL_div_covs(cov_second, cov_first)))
