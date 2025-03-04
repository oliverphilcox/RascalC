"This reads cosmodesi/pycorr .npy file(s) and generates input xi text file for RascalC to use"

import sys

## PARAMETERS
if len(sys.argv) < 5:
    print("Usage: python convert_xi_from_pycorr.py {INPUT_NPY_FILE1} [{INPUT_NPY_FILE2} ...] {OUTPUT_XI_DAT_FILE} {R_STEP} {N_MU}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.pycorr_utils.input_xi import convert_xi_from_pycorr_files

infile_names = sys.argv[1:-3]
outfile_name = str(sys.argv[-3])
r_step = float(sys.argv[-2])
n_mu = int(sys.argv[-1])

convert_xi_from_pycorr_files(infile_names, outfile_name, n_mu, r_step)