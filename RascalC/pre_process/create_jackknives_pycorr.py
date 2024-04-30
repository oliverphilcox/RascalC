""" This function will assign a jackknife region to each input particle, in a way compatible with pycorr used in DESI data processing. Data is saved as a 5 column text file."""

import numpy as np
import os
import pycorr


def get_subsampler_xirunpc(ref_positions: np.ndarray[float], njack: int, position_type: str = 'rdd') -> pycorr.KMeansSubsampler:
    if len(ref_positions) != 3: ref_positions = ref_positions.T
    if len(ref_positions) != 3: raise ValueError("The reference positions do not appear 3-dimensional")
    # Remove the OMP_ environment variables temporarily
    omp_keys = [key for key in os.environ if key.startswith("OMP_")]
    omp_backup = {key: os.environ.pop(key) for key in omp_keys}
    os.environ["OMP_NUM_THREADS"] = "1"
    # Create the subsampler. The above steps disable multi-threading, which is known to interfere with RascalC C++ code. Later usage of the subsampler should remain single-threaded as we want it. The jackknife assignment takes very little time compared to the pair counting and especially RascalC importance sampling.
    subsampler = pycorr.KMeansSubsampler('angular', positions = ref_positions, position_type = position_type, nsamples = njack, nside = 512, random_state = 42)
    # Restore the OMP_ environment variables modified above
    os.environ.pop("OMP_NUM_THREADS")
    for key in omp_keys: os.environ[key] = omp_backup[key]
    return subsampler


def add_subsample_labels(subsampler: pycorr.twopoint_jackknife.BaseSubsampler, pos_etc: np.ndarray[float], position_type: str = "xyz") -> np.ndarray[float]:
    positions = pos_etc[:, :3]
    labels = subsampler.label(positions.T, position_type = position_type)
    return np.column_stack((pos_etc, labels))


def create_jackknives_pycorr_files(reffile_name: str, infile_name: str, outfile_name: str, njack: int, print_function = print):
    import time
    init = time.time()

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

    end = time.time() - init
    print_function('Task completed in %.2f seconds' %end)
