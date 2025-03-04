## Script to perform an extra convergence check on full integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Should work in any case - default, jackknife, Legendre, multi-tracer - as it utilizes universal data from RascalC file

import sys

# PARAMETERS
if len(sys.argv) not in (2, 3):
    print("Usage: python convergence_check_extra.py {RASCALC_RESULTS_FILE} [{N_SUBSAMPLES_TO_USE}]")
    sys.exit(1)

from utils import adjust_path, get_arg_safe
adjust_path()
from RascalC.convergence_check_extra import convergence_check_extra_file

rascalc_results = str(sys.argv[1])
n_samples = get_arg_safe(2, int)

convergence_check_extra_file(rascalc_results, n_samples, print_function = print)