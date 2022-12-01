## Script to perform an extra convergence check on full integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Should work in any case - default, jackknife, Legendre, multi-tracer - as it utilizes universal data from RascalC file

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv) not in (2, 3):
    print("Usage: python convergence_check_extra.py {RASCALC_RESULTS_FILE} [{N_SUBSAMPLES_TO_USE}]")
    sys.exit()

rascalc_results = str(sys.argv[1])

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

with np.load(rascalc_results) as f:
    c_samples = f["individual_theory_covariances"]
    jack_key = "individual_theory_jackknife_covariances"
    do_jack = (jack_key in f.keys())
    if do_jack: cj_samples = f[jack_key]

n_samples = int(sys.argv[2]) if len(sys.argv) >= 3 else len(c_samples)
n_samples_2 = n_samples // 2

def cmp_cov_print(cov_first, cov_second):
    print("RMS eigenvalues of inverse tests for cov half-estimates are %.2e and %.2e" % (rms_eig_inv_test_covs(cov_first, cov_second), rms_eig_inv_test_covs(cov_second, cov_first)))
    print("KL divergences between cov half-estimates are %.2e and %.2e" % (KL_div_covs(cov_first, cov_second), KL_div_covs(cov_second, cov_first)))

def cmp_cov_splittings(c_samples):
    print("First splitting")
    cov_first = np.mean(c_samples[:n_samples_2], axis=0)
    cov_second = np.mean(c_samples[n_samples_2:n_samples], axis=0)
    cmp_cov_print(cov_first, cov_second)

    print("Second splitting")
    cov_first = np.mean(c_samples[:n_samples:2], axis=0)
    cov_second = np.mean(c_samples[1:n_samples:2], axis=0)
    cmp_cov_print(cov_first, cov_second)

print("Full covariance")
cmp_cov_splittings(c_samples)
if do_jack:
    print("Jack covariance")
    cmp_cov_splittings(cj_samples)