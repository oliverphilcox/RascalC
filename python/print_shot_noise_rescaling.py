## Simple convenience script to print shot noise rescaling from one or many RascalC files

import sys
import numpy as np

# PARAMETERS
if len(sys.argv) < 2:
    print("Usage: python print_shot_noise_rescaling.py {RASCALC_RESULTS_1} [{RASCALC_RESULTS_2} ...]")
    sys.exit(1)

for filename in sys.argv[1:]:
    with np.load(filename) as f:
        print(f"{filename}: {f['shot_noise_rescaling'] = }")