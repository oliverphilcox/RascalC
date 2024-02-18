"This reads cosmodesi/pycorr .npy file(s) and generates sample covariance of xi_l(s) in text format"

import pycorr
import numpy as np
from .utils import reshape_pycorr


def sample_cov_multipoles_from_pycorr(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], max_l: int, r_step: float = 1, r_max: float = np.inf):
    if max_l < 0: raise ValueError("Maximal multipole can not be negative")
    assert max_l % 2 == 0, "Only even multipoles supported"
    if len(xi_estimators) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(xi_estimators[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(xi_estimators_c) != len(xi_estimators[0]) for xi_estimators_c in xi_estimators[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    ells = np.arange(0, max_l+1, 2)
    # convert each xi estimator to multipoles array, and then turn list of lists into array too
    xi = np.array([[reshape_pycorr(xi_estimator, r_step = r_step, r_max = r_max, n_mu = None).get_corr(mode='poles', ells = ells) for xi_estimator in xi_estimators_c] for xi_estimators_c in xi_estimators])
    # now indices are [c, s, l, r]: correlation number, sample number, multipole index and then radial bin
    # need [s, c, l, r]
    xi = xi.transpose(1, 0, 2, 3)
    # now flatten all dimensions except the samples
    xi = xi.reshape(xi.shape[0], -1)
    return np.cov(xi.T) # xi has to be transposed, because variables (bins) are in columns (2nd index) of it and np.cov expects otherwise.
    # Weights are assumed the same, hard to figure out alternatives, and they do not seem necessary


def sample_cov_multipoles_from_pycorr_to_file(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], outfile_name, max_l: int, r_step: float = 1, r_max: float = np.inf):
    np.savetxt(outfile_name, sample_cov_multipoles_from_pycorr(xi_estimators, max_l, r_step, r_max))


def sample_cov_multipoles_from_pycorr_files(infile_names: list[list[str]], outfile_name: str, max_l: int, r_step: float = 1, r_max: float = np.inf):
    if max_l < 0: raise ValueError("Maximal multipole can not be negative")
    assert max_l % 2 == 0, "Only even multipoles supported"
    if len(infile_names) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(infile_names[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(infile_names_c) != len(infile_names[0]) for infile_names_c in infile_names[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    xi_estimators = [[pycorr.TwoPointCorrelationFunction.load(infile_name) for infile_name in infile_names_c] for infile_names_c in infile_names]
    sample_cov_multipoles_from_pycorr_to_file(xi_estimators, outfile_name, max_l, r_step, r_max)
