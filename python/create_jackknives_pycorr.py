""" This function will assign a jackknife region to each input particle, in a way compatible with pycorr used in DESI data processing. Data is saved as a 5 column text file."""

import numpy as np
import pycorr
import sys


def get_subsampler_xirunpc(ref_positions: np.ndarray[float], njack: int) -> pycorr.KMeansSubsampler:
    if len(ref_positions) != 3: ref_positions = ref_positions.T
    if len(ref_positions) != 3: raise ValueError("The reference positions do not appear 3-dimensional")
    return pycorr.KMeansSubsampler('angular', positions = ref_positions, position_type = 'rdd', nsamples = njack, nside = 512, random_state = 42)

def add_subsample_labels(subsampler: pycorr.twopoint_jackknife.BaseSubsampler, pos_etc: np.ndarray[float], position_type: str = "xyz") -> np.ndarray[float]:
    positions = pos_etc[:, :3]
    labels = subsampler.label(positions.T, position_type = position_type)
    return np.column_stack((pos_etc, labels))

def create_jackknives_pycorr_files(reffile_name: str, infile_name: str, outfile_name: str, njack: int, print_function = print):
    import time
    init=time.time()

    print_function("Reading reference positions from %s" % reffile_name)
    ref_positions = np.loadtxt(reffile_name, usecols=range(3)) # only read positions

    print_function("Initializing the subsampler")
    subsampler = get_subsampler_xirunpc(ref_positions, njack)

    print_function("Reading positions and weights from %s" % infile_name)
    pos_etc = np.loadtxt(infile_name)

    print_function("Assigning jackknives")
    pos_weights_jack = add_subsample_labels(subsampler, pos_etc)

    print_function("Saving to %s" % outfile_name)
    np.savetxt(outfile_name, pos_weights_jack)

    end=time.time()-init
    print_function('Task completed in %.2f seconds' %end)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) != 5:
        print("Usage: python create_jackknives_pycorr.py {REF_RDZ_FILE} {INPUT_XYZW_FILE} {OUTPUT_XYZWJ_FILE} {NJACK}.")
        sys.exit(1)
    reffile_name = str(sys.argv[1])
    infile_name = str(sys.argv[2])
    outfile_name = str(sys.argv[3])
    njack = int(sys.argv[4])

    create_jackknives_pycorr_files(reffile_name, infile_name, outfile_name, njack)