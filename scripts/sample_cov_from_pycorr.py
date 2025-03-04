"This reads cosmodesi/pycorr .npy file(s) and generates sample covariance of xi(s,mu) in text format"

import sys

## PARAMETERS
if len(sys.argv) < 8:
    print("Usage: python sample_cov_from_pycorr.py {INPUT_NPY_FILE1} {INPUT_NPY_FILE2} [{INPUT_NPY_FILE3} ...] {OUTPUT_COV_FILE} {R_STEP} {N_MU} {R_MAX} {N_CORRELATIONS}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.pycorr_utils.sample_cov import sample_cov_from_pycorr_files

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