"This reads a cosmodesi/pycorr .npy file and generates binned pair counts text file for RascalC to use"

import pycorr
import numpy as np
from .utils import reshape_pycorr


def get_counts_from_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator, counts_factor: float | None = None, split_above: float = np.inf) -> np.typing.NDArray[np.float64]:
    if not counts_factor: # use normalized counts
        return xi_estimator.R1R2.normalized_wcounts()
    paircounts = xi_estimator.R1R2.wcounts / counts_factor
    nonsplit_mask = (xi_estimator.sepavg(axis=0) < split_above)
    if split_above > 0: paircounts[nonsplit_mask] /= counts_factor # divide once more below the splitting scale
    return paircounts


def convert_counts_from_pycorr_to_file(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator, outfile_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf):
    xi_estimator_reshaped = reshape_pycorr(xi_estimator, n_mu, r_step, r_max)
    paircounts = get_counts_from_pycorr(xi_estimator_reshaped, counts_factor, split_above)
    ## Write to file using numpy funs
    np.savetxt(outfile_name, paircounts.reshape(-1, 1)) # the file always has 1 column


def convert_counts_from_pycorr_files(infile_name: str, outfile_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf):
    xi_estimator_orig = pycorr.TwoPointCorrelationFunction.load(infile_name)
    convert_counts_from_pycorr_to_file(xi_estimator_orig, outfile_name, n_mu, r_step, r_max, counts_factor, split_above)
