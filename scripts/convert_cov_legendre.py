"Slightly more complicated (than convert_cov.py) convenience script that reads RascalC Legendre mode results and saves full cov to text file, changing the indexing from [r, l] to [l, r], in addition checking eigenvalues of bias matrix"

import sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python convert_cov_legendre.py {RASCALC_RESULTS_FILE} {MAX_L} {OUTPUT_COV_FILE}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.cov_utils import export_cov_legendre

rascalc_results = str(sys.argv[1])
max_l = int(sys.argv[2])
output_cov_file = str(sys.argv[3])

export_cov_legendre(rascalc_results, max_l, output_cov_file)
