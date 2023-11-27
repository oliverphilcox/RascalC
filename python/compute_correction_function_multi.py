### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a two surveys. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import sys
import os
import numpy as np
from .compute_correction_function import compute_V_n_w_bar_from_file, compute_phi_periodic, load_RR, compute_phi_aperiodic

def compute_correction_function_multi(random_file: str, random_file2: str, binfile: str, outdir: str, periodic: bool, RR_file: str | None = None, RR_file12: str | None = None, RR_file2: str | None = None, print_function = print) -> None:
    if periodic:
        print_function("Assuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere")
    elif any(file is None for file in (RR_file, RR_file12, RR_file2)):
        raise TypeError("All the RR files must be specified if aperiodic")

    V1, n_bar1, w_bar1 = compute_V_n_w_bar_from_file(random_file, print_function = print_function)
    V2, n_bar2, w_bar2 = compute_V_n_w_bar_from_file(random_file2, index = 2, print_function = print_function)

    # Load in binning files
    r_bins = np.loadtxt(binfile)
    n = len(r_bins)

    if periodic:
        ## Periodic correction function is simple

        phi_11 = compute_phi_periodic(w_bar1, w_bar1, n_bar1, n_bar1, V1, n)
        phi_12 = compute_phi_periodic(w_bar1, w_bar2, n_bar1, n_bar2, np.sqrt(V1*V2), n)
        phi_22 = compute_phi_periodic(w_bar2, w_bar2, n_bar2, n_bar2, V2, n)
    else:
        ## Continue for aperiodic case

        # Load RR counts
        RR = load_RR(RR_file, n)
        m = RR.shape[1]

        RR2 = load_RR(RR_file2, n)
        if RR2.shape != RR.shape: raise ValueError("Need the same binning for all RR files")

        RR12 = load_RR(RR_file12, n)
        if RR12.shape != RR.shape: raise ValueError("Need the same binning for all RR files")

        # Now run this
        phi_11 = compute_phi_aperiodic(w_bar1, w_bar1, n_bar1, n_bar1, V1, r_bins, RR, "11", print_function)
        phi_12 = compute_phi_aperiodic(w_bar1, w_bar2, n_bar1, n_bar2, np.sqrt(V1*V2), r_bins, RR, "12", print_function) # use geometrically average volume, should be not critical
        phi_12 = compute_phi_aperiodic(w_bar2, w_bar2, n_bar2, n_bar2, V2, RR, r_bins, "22", print_function)

    roots = ['11', '12', '22']
    phis = [phi_11, phi_12, phi_22]

    for phi, index in zip(phis, roots):
        if periodic:
            outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_periodic_%s.txt'%(n, index))
        else:
            outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_m%d_11.%s'%(n, m, index))
        np.savetxt(outfile, phi)
        print_function("Saved (normalized) output for field %s to %s"%(index, outfile))

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (6, 9):
        print("Usage: python compute_correction_function_multi.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS_11} {RR_COUNTS_12} {RR_COUNTS_22}]")
        sys.exit(1)
    random_file = str(sys.argv[1])
    random_file2 = str(sys.argv[2])
    binfile = str(sys.argv[3])
    outdir = str(sys.argv[4])
    periodic = int(sys.argv[5])

    from .utils import get_arg_safe
    RR_file = get_arg_safe(6)
    RR_file12 = get_arg_safe(7)
    RR_file2 = get_arg_safe(8)

    compute_correction_function_multi(random_file, random_file2, binfile, outdir, periodic, RR_file, RR_file12, RR_file2)