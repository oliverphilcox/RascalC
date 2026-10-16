"This reads cosmodesi/pycorr .npy file(s) and generates input xi text file for RascalC to use"

import argparse

parser = argparse.ArgumentParser(description="This script reads cosmodesi/pycorr .npy file(s) and generates input xi text file for RascalC to use")
parser.add_argument("pycorr_files", type=str, nargs='+', help="pycorr .npy filename(s). if multiple are provided, the correlation function will be averaged")
parser.add_argument("output_xi_file", type=str, help="output correlation function (xi) text file")
parser.add_argument("r_step", type=float, help="desired width of radial bins")
parser.add_argument("n_mu_bins", type=int, help="number of angular (mu) bins")
args = parser.parse_args()

from utils import adjust_path
adjust_path()
from RascalC.pycorr_utils.input_xi import convert_xi_from_pycorr_files

convert_xi_from_pycorr_files(args.pycorr_files, args.outfile_name, args.n_mu_bins, args.r_step)