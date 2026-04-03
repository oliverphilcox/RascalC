"This generates jackknife xi, weight, paircounts, and total paircounts text files for RascalC from pycorr TwoPointCorrelationFunction objects/files"

import pycorr
import numpy as np
import numpy.typing as npt
from .utils import reshape_pycorr, fix_bad_bins_pycorr, write_xi_file
from .counts import get_counts_from_pycorr


def jack_realization_rascalc(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, i: int) -> pycorr.TwoPointEstimator:
    # returns RascalC-framed jackknife realization, different from implemented in pycorr
    cls = jack_estimator.__class__.__bases__[0]
    kw = {}
    for name in jack_estimator.count_names:
        counts = getattr(jack_estimator, name)
        kw[name] = counts.auto[i] + 0.5 * (counts.cross12[i] + counts.cross21[i]) # j1 x all2 + all1 x j2 = 2 x j1 x j2 + j1 x (all2 - j2) + j2 x (all1 - j1); by conventions from jackknife_weights{,_cross}.py needs to be divided by 2
    return fix_bad_bins_pycorr(cls(**kw))


def jack_realizations_rascalc(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator) -> list[pycorr.TwoPointEstimator]:
    return [jack_realization_rascalc(jack_estimator, i) for i in jack_estimator.realizations]


def get_jack_xi_weights_counts_from_pycorr(jack_estimator: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, counts_factor: float | None = None, split_above: float = np.inf) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    realizations = jack_realizations_rascalc(jack_estimator)

    xi_jack = np.array([jack.corr.ravel() for jack in realizations]) # already wrapped
    jack_pairs = np.array([get_counts_from_pycorr(jack, counts_factor, split_above).ravel() for jack in realizations]) # already wrapped
    jack_weights = jack_pairs / np.sum(jack_pairs, axis=0)[None, :] # weights are pair counts normalized by their total
    return xi_jack, jack_weights, jack_pairs


def convert_jack_xi_weights_counts_from_pycorr_to_files(xi_estimator_orig: pycorr.twopoint_jackknife.JackknifeTwoPointEstimator, xi_jack_name: str, jackweights_name: str, jackpairs_name: str, binpairs_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf) -> None:
    xi_estimator = reshape_pycorr(xi_estimator_orig, n_mu, r_step, r_max)
    binpairs = get_counts_from_pycorr(xi_estimator, counts_factor, split_above).ravel()

    xi_jack, jack_weights, jack_pairs = get_jack_xi_weights_counts_from_pycorr(xi_estimator, counts_factor, split_above)
    if not np.allclose(np.sum(jack_pairs, axis=0), binpairs): raise ValueError("Total counts mismatch")

    ## Write to files using numpy functions
    np.savetxt(binpairs_name, binpairs.reshape(-1, 1)) # this file must have 1 column
    write_xi_file(xi_jack_name, xi_estimator.sepavg(axis = 0), xi_estimator.sepavg(axis = 1), xi_jack)
    jack_numbers = np.array(xi_estimator.realizations).reshape(-1, 1) # column of jackknife numbers, may be useless but needed for format compatibility
    np.savetxt(jackweights_name, np.column_stack((jack_numbers, jack_weights)))
    np.savetxt(jackpairs_name, np.column_stack((jack_numbers, jack_pairs)))


def convert_jack_xi_weights_counts_from_pycorr_files(infile_name: str, xi_jack_name: str, jackweights_name: str, jackpairs_name: str, binpairs_name: str, n_mu: int | None = None, r_step: float = 1, r_max: float = np.inf, counts_factor: float | None = None, split_above: float = np.inf) -> None:
    xi_estimator = pycorr.TwoPointCorrelationFunction.load(infile_name)
    convert_jack_xi_weights_counts_from_pycorr_to_files(xi_estimator, xi_jack_name, jackweights_name, jackpairs_name, binpairs_name, n_mu, r_step, r_max, counts_factor, split_above)
