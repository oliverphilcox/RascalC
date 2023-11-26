"This reads cosmodesi/pycorr .npy file(s) and generates sample covariance of xi(s,mu) in text format"

import pycorr
import numpy as np
import sys
from utils import reshape_pycorr


def sample_cov_from_pycorr(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf):
    if len(xi_estimators) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(xi_estimators[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(xi_estimators_c) != len(xi_estimators[0]) for xi_estimators_c in xi_estimators[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    # convert each xi estimator to binned xi array, and then turn list of lists into array too
    xi = np.array([[reshape_pycorr(xi_estimator, n_mu, r_step, r_max).corr for xi_estimator in xi_estimators_c] for xi_estimators_c in xi_estimators])
    # now indices are [c, s, r, m]: correlation number, sample number, radial bin and then mu bin
    # need [s, c, r, m]
    xi = xi.transpose(1, 0, 2, 3)
    # now flatten all dimensions except the samples
    xi = xi.reshape(xi.shape[0], -1)
    return np.cov(xi.T) # xi has to be transposed, because variables (bins) are in columns (2nd index) of it and np.cov expects otherwise.
    # Weights are assumed the same, hard to figure out alternatives, and they do not seem necessary

def sample_cov_from_pycorr_to_file(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], outfile_name, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf):
    np.savetxt(outfile_name, sample_cov_from_pycorr(xi_estimators, n_mu, r_step, r_max))

def sample_cov_from_pycorr_files(infile_names: list[list[str]], outfile_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf):
    if len(infile_names) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(infile_names[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(infile_names_c) != len(infile_names[0]) for infile_names_c in infile_names[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    xi_estimators = [[pycorr.TwoPointCorrelationFunction.load(infile_name) for infile_name in infile_names_c] for infile_names_c in infile_names]
    sample_cov_from_pycorr_to_file(xi_estimators, outfile_name, n_mu, r_step, r_max)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) < 8:
        print("Usage: python sample_cov_from_pycorr.py {INPUT_NPY_FILE1} {INPUT_NPY_FILE2} [{INPUT_NPY_FILE3} ...] {OUTPUT_COV_FILE} {R_STEP} {N_MU} {R_MAX} {N_CORRELATIONS}.")
        sys.exit(1)
    infile_names = sys.argv[1:-5]
    outfile_name = str(sys.argv[-5])
    r_step = float(sys.argv[-4])
    n_mu = int(sys.argv[-3])
    r_max = float(sys.argv[-2])
    n_corr = int(sys.argv[-1])

    assert n_corr >= 1, "Need to have at least one correlation"
    n_files = len(infile_names)
    assert n_files % n_corr == 0, "Need to have the same number of files for all correlations"
    n_samples = n_files // n_corr
    infile_names = [infile_names[n_samples * i, n_samples * (i+1)] for i in range(n_corr)]

    sample_cov_from_pycorr_files(infile_names, outfile_name, n_mu, r_step, r_max)