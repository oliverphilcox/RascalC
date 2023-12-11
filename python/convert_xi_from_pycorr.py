"This reads cosmodesi/pycorr .npy file(s) and generates input xi text file for RascalC to use"

import pycorr
import numpy as np
import sys
from .utils import reshape_pycorr, write_xi_file


def get_input_xi_from_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator):
    # assume already wrapped; for input xi need to divide by SS instead of RR in post-recon case, in pre-recon case RR=SS so should work too
    return xi_estimator.corr * xi_estimator.R1R2.normalized_wcounts() / xi_estimator.S1S2.normalized_wcounts()

def convert_xi_from_pycorr_to_file(xi_estimators: list[pycorr.twopoint_estimator.BaseTwoPointEstimator], outfile_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, print_function = print) -> tuple[float, float]:
    # Compute mean data sizes
    mean_data_size1 = np.mean([xi_estimator.D1D2.size1 for xi_estimator in xi_estimators])
    mean_data_size2 = np.mean([xi_estimator.D1D2.size2 for xi_estimator in xi_estimators])

    print_function(f"Mean size of data 1 is {mean_data_size1:.6e}")
    print_function(f"Mean size of data 2 is {mean_data_size2:.6e}")
    np.savetxt(outfile_name + ".ndata", np.array((mean_data_size1, mean_data_size2))) # save them for later

    # Reshape the estimators (includes fixing)
    xi_estimators = [reshape_pycorr(xi_estimator, n_mu, r_step, r_max) for xi_estimator in xi_estimators]

    # Sum the estimators to get total
    xi_estimator = sum(xi_estimators)

    ## Write to file using numpy funs
    write_xi_file(outfile_name, xi_estimator.sepavg(axis=0), xi_estimator.sepavg(axis=1), get_input_xi_from_pycorr(xi_estimator))

    return mean_data_size1, mean_data_size2

def convert_xi_from_pycorr_files(infile_names: list[str], outfile_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, print_function = print):
    # load all pycorr files
    xi_estimators = [pycorr.TwoPointCorrelationFunction.load(infile_name) for infile_name in infile_names]
    return convert_xi_from_pycorr_to_file(xi_estimators, outfile_name, n_mu, r_step, r_max, print_function)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) < 5:
        print("Usage: python convert_xi_from_pycorr.py {INPUT_NPY_FILE1} [{INPUT_NPY_FILE2} ...] {OUTPUT_XI_DAT_FILE} {R_STEP} {N_MU}.")
        sys.exit(1)
    infile_names = sys.argv[1:-3]
    outfile_name = str(sys.argv[-3])
    r_step = float(sys.argv[-2])
    n_mu = int(sys.argv[-1])

    convert_xi_from_pycorr_files(infile_names, outfile_name, n_mu, r_step)