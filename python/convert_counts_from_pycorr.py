"This reads a cosmodesi/pycorr .npy file and generates binned pair counts text file for RascalC to use"

import pycorr
import numpy as np
import sys
from utils import reshape_pycorr


def get_counts_from_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator, counts_factor: float | None = None, split_above: float = np.inf) -> np.ndarray[float]:
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

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) not in (3, 4, 5, 6, 7, 8):
        print("Usage: python convert_counts_from_pycorr.py {INPUT_NPY_FILE} {OUTPUT_PAIRCOUNTS_TEXT_FILE} [{R_STEP} [{N_MU} [{COUNTS_FACTOR} [{SPLIT_ABOVE}] [{R_MAX}]]]]]. COUNTS_FACTOR=0 is a special value to use normalized counts.")
        sys.exit(1)
    infile_name = str(sys.argv[1])
    outfile_name = str(sys.argv[2])
    from utils import get_arg_safe
    r_step = get_arg_safe(3, float, 1)
    n_mu = get_arg_safe(4, int, None)
    counts_factor = get_arg_safe(5, float, None) # basically number of randoms used for these counts, used to convert from total to 1 catalog count estimate
    split_above = get_arg_safe(6, float, np.inf) # divide weighted RR counts by counts_factor**2 below this and by counts_factor above
    r_max = get_arg_safe(7, float, np.inf) # if given, limit used r_bins at that

    if n_mu == 0: n_mu = None
    if counts_factor == 0: counts_factor = None
    convert_counts_from_pycorr_files(infile_name, outfile_name, n_mu, r_step, r_max, counts_factor, split_above)