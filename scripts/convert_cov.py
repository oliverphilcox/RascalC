"Simple convenience script that reads RascalC results and saves full cov to text file, in addition checking eigenvalues of bias matrix"

import sys

## PARAMETERS
if len(sys.argv) != 3:
    print("Usage: python convert_cov.py {RASCALC_RESULTS_FILE} {OUTPUT_COV_FILE}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.cov_utils import export_cov

rascalc_results = str(sys.argv[1])
output_cov_file = str(sys.argv[2])

export_cov(rascalc_results, output_cov_file)