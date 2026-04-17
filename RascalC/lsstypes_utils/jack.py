"This generates jackknife xi, weight, paircounts, and total paircounts text files for RascalC from lsstypes Count2JackknifeCorrelation objects/files"

import lsstypes
import numpy as np
import numpy.typing as npt
from .counts import get_counts_from_lsstypes


def jack_realization_rascalc(jack_estimator: lsstypes.Count2JackknifeCorrelation, ii: int) -> lsstypes.Count2Correlation:
    # returns RascalC-framed jackknife realization, different from implemented in lsstypes
    # based on lsstypes.Count2JackknifeCorrelation.realization and lsstypes.Count2Jackknife.realization
    kw = {}
    for count_name in jack_estimator.count_names:
        jack_counts : lsstypes.Count2Jackknife = jack_estimator.get(count_name)
        counts = jack_counts.get(realization=ii, cross='ii').values('counts') + 0.5 * (jack_counts.get(realization=ii, cross='ij').values('counts') + jack_counts.get(realization=ii, cross='ji').values('counts')) # j1 x all2 + all1 x j2 = 2 x j1 x j2 + j1 x (all2 - j2) + j2 x (all1 - j1); by conventions from jackknife_weights{,_cross}.py needs to be divided by 2
        kw[count_name] = lsstypes.Count2(counts=counts, norm=jack_counts.values('norm'), attrs=jack_counts.get(realization=ii, cross='ii').attrs) # we want the norm to be the same as for the total counts for further RascalC usage, otherwise the normalized counts of the jackknife realizations would not sum to the total normalized counts
        # could set the proper coordinates (average separations) as in lsstypes.Count2Jackknife.realization, but RascalC is not going to use them
        # could also set size attributes as in lsstypes.Count2Jackknife.realization, but it is not clear how to do this right, and RascalC is not going to use them either
    return lsstypes.Count2Correlation(**kw, estimator=jack_estimator.estimator, attrs=jack_estimator.attrs)


def jack_realizations_rascalc(jack_estimator: lsstypes.Count2JackknifeCorrelation) -> list[lsstypes.Count2Correlation]:
    return [jack_realization_rascalc(jack_estimator, ii) for ii in jack_estimator.realizations]


def get_jack_xi_weights_counts_from_lsstypes(jack_estimator: lsstypes.Count2JackknifeCorrelation, counts_factor: float | None = None, split_above: float = np.inf) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    realizations = jack_realizations_rascalc(jack_estimator)

    xi_jack = np.array([jack.value().ravel() for jack in realizations]) # already wrapped
    jack_pairs = np.array([get_counts_from_lsstypes(jack, counts_factor, split_above).ravel() for jack in realizations]) # already wrapped
    jack_weights = jack_pairs / np.sum(jack_pairs, axis=0)[None, :] # weights are pair counts normalized by their total
    return xi_jack, jack_weights, jack_pairs
