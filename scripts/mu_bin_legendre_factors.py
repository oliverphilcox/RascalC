"Simple script to produce the file of mu bin Legendre factors for the C++ code in LEGENDRE_MIX mode."

import argparse

parser = argparse.ArgumentParser(description="Simple script to produce the file of mu bin Legendre factors for the C++ code in LEGENDRE_MIX mode.")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
parser.add_argument("max_l", type=int, help="maximum multipole index")
parser.add_argument("output_dir", type=str, help="directory to save the resulting file")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.mu_bin_legendre_factors import write_mu_bin_legendre_factors

write_mu_bin_legendre_factors(args.n_mu_bins, args.max_l, args.output_dir)