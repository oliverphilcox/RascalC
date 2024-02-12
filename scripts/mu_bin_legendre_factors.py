"Simple script to produce the file of mu bin Legendre factors for the C++ code in LEGENDRE_MIX mode."

import sys

## PARAMETERS
if len(sys.argv) != 4:
    print("Usage: python mu_bin_legendre_factors.py {N_MU_BINS} {MAX_L} {OUTPUT_DIR}.")
    sys.exit(1)

from utils import adjust_path
adjust_path()
from RascalC.mu_bin_legendre_factors import write_mu_bin_legendre_factors

n_mu_bins = int(sys.argv[1])
max_l = int(sys.argv[2])
output_dir = str(sys.argv[3])

write_mu_bin_legendre_factors(n_mu_bins, max_l, output_dir)