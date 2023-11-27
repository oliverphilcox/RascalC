"This reads a cosmodesi/pycorr .npy file and generates jackknife xi, weight, paircounts, and total paircounts text files for RascalC to use"

import pycorr
import numpy as np
import sys
from .utils import reshape_pycorr, fix_bad_bins_pycorr, write_xi_file
from .convert_counts_from_pycorr import get_counts_from_pycorr


def jack_realization_rascalc(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, i) -> pycorr.TwoPointEstimator:
    # returns RascalC-framed jackknife realization, different from implemented in pycorr
    cls = jack_estimator.__class__.__bases__[0]
    kw = {}
    for name in jack_estimator.count_names:
        counts = getattr(jack_estimator, name)
        kw[name] = counts.auto[i] + 0.5 * (counts.cross12[i] + counts.cross21[i]) # j1 x all2 + all1 x j2 = 2 x j1 x j2 + j1 x (all2 - j2) + j2 x (all1 - j1); by conventions from jackknife_weights{,_cross}.py needs to be divided by 2
    return fix_bad_bins_pycorr(cls(**kw))

def jack_realizations_rascalc(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator) -> list[pycorr.TwoPointEstimator]:
    return [jack_realization_rascalc(jack_estimator, i) for i in jack_estimator.realizations]

def get_jack_xi_weights_counts_from_pycorr(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, counts_factor: float | None = None, split_above: float = np.inf) -> tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]]:
    realizations = jack_realizations_rascalc(jack_estimator)

    xi_jack = np.array([jack.corr.ravel() for jack in realizations]) # already wrapped
    if counts_factor: # nonzero value
        jack_pairs = np.array([jack.R1R2.wcounts.ravel() for jack in realizations]) / counts_factor # already wrapped
        nonsplit_mask = (jack_estimator.sepavg(axis=0) < split_above)
        if split_above > 0: jack_pairs[:, nonsplit_mask] /= counts_factor # divide once more below the splitting scale
    else: # zero value, use normalized counts
        jack_pairs = np.array([(jack.R1R2.wcounts / jack_estimator.R1R2.wnorm).ravel() for jack in realizations]) # already wrapped
    jack_weights = jack_pairs / np.sum(jack_pairs, axis=0)[None, :] # weights are pair counts normalized by their total
    return xi_jack, jack_weights, jack_pairs

def convert_jack_xi_weights_counts_from_pycorr_to_files(xi_estimator_orig: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, xi_jack_name: str, jackweights_name: str, jackpairs_name: str, binpairs_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf):
    xi_estimator = reshape_pycorr(xi_estimator_orig, n_mu, r_step, r_max)
    binpairs = get_counts_from_pycorr(xi_estimator, counts_factor, split_above)

    xi_jack, jack_weights, jack_pairs = get_jack_xi_weights_counts_from_pycorr(xi_estimator, counts_factor, split_above)
    if not np.allclose(np.sum(jack_pairs, axis=0), binpairs): raise ValueError("Total counts mismatch")

    ## Write to files using numpy functions
    np.savetxt(binpairs_name, binpairs.reshape(-1, 1)) # this file must have 1 column
    write_xi_file(xi_jack_name, xi_estimator.sepavg(axis = 0), xi_estimator.sepavg(axis = 1), xi_jack)
    jack_numbers = np.array(xi_estimator.realizations).reshape(-1, 1) # column of jackknife numbers, may be useless but needed for format compatibility
    np.savetxt(jackweights_name, np.hstack((jack_numbers, jack_weights)))
    np.savetxt(jackpairs_name, np.hstack((jack_numbers, jack_pairs)))

def convert_jack_xi_weights_counts_from_pycorr_files(infile_name: str, xi_jack_name: str, jackweights_name: str, jackpairs_name: str, binpairs_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf):
    xi_estimator = pycorr.TwoPointCorrelationFunction.load(infile_name)
    convert_jack_xi_weights_counts_from_pycorr_to_files(xi_estimator, xi_jack_name, jackweights_name, jackpairs_name, binpairs_name, n_mu, r_step, r_max, counts_factor, split_above)

if __name__ == "__main__": # if invoked as a script
    ## PARAMETERS
    if len(sys.argv) not in (6, 7, 8, 9, 10, 11):
        print("Usage: python convert_xi_jack_from_pycorr.py {INPUT_NPY_FILE} {OUTPUT_XI_JACK_FILE} {OUTPUT_JACKKNIFE_WEIGHTS_FILE} {OUTPUT_JACKKNIFE_PAIRCOUNTS_FILE} {OUTPUT_BINNED_PAIRCOUNTS_FILE} [{R_STEP} [{N_MU} [{COUNTS_FACTOR} [{SPLIT_ABOVE} [{R_MAX}]]]]]. COUNTS_FACTOR=0 is a special value to use normalized counts.")
        sys.exit(1)
    infile_name = str(sys.argv[1])
    xi_jack_name = str(sys.argv[2])
    jackweights_name = str(sys.argv[3])
    jackpairs_name = str(sys.argv[4])
    binpairs_name = str(sys.argv[5])
    from .utils import get_arg_safe
    r_step = get_arg_safe(6, float, 1)
    n_mu = get_arg_safe(7, int, None)
    counts_factor = get_arg_safe(8, float, None) # basically number of randoms used for these counts, used to convert from total to 1 catalog count estimate
    split_above = get_arg_safe(9, float, np.inf) # divide weighted RR counts by counts_factor**2 below this and by counts_factor above
    r_max = get_arg_safe(10, float, np.inf) # if given, limit used r_bins at that

    convert_jack_xi_weights_counts_from_pycorr_files(infile_name, xi_jack_name, jackweights_name, jackpairs_name, binpairs_name, n_mu, r_step, r_max, counts_factor, split_above)