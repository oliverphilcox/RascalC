"Universal functions for pycorr TwoPointEstimator and lsstypes Count2Correlation, to avoid code duplication in interface.py and other scripts."
import numpy as np
import numpy.typing as npt
import pycorr
import lsstypes
from typing import Callable, Any
from .pycorr_utils.utils import fix_bad_bins_pycorr
from .lsstypes_utils.utils import get_edges_from_lsstypes_to_pycorr, fix_bad_bins_lsstypes, wrap_correlation


def allcount_switch_function(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation, func_pycorr: Callable[[pycorr.twopoint_estimator.BaseTwoPointEstimator], Any], func_lsstypes: Callable[[lsstypes.Count2Correlation], Any]) -> Any:
    "switch function to apply func_pycorr if allcounts is a pycorr estimator and func_lsstypes if allcounts is a lsstypes Count2Correlation"
    if isinstance(allcounts, pycorr.twopoint_estimator.BaseTwoPointEstimator):
        return func_pycorr(allcounts)
    if isinstance(allcounts, lsstypes.Count2Correlation):
        return func_lsstypes(allcounts)
    raise TypeError("allcounts must be either a pycorr BaseTwoPointEstimator or a lsstypes Count2Correlation")


def get_s_edges_from_allcounts(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation) -> npt.NDArray[np.float64]:
    "get separation/radial bin edges from allcounts or xi estimator"
    return allcount_switch_function(allcounts, lambda x: x.edges[0], lambda x: get_edges_from_lsstypes_to_pycorr(x, 's'))


def get_mu_edges_from_allcounts(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation) -> npt.NDArray[np.float64]:
    "get mu/angular bin edges from allcounts or xi estimator"
    return allcount_switch_function(allcounts, lambda x: x.edges[1], lambda x: get_edges_from_lsstypes_to_pycorr(x, 'mu'))


def get_s_avg_from_allcounts(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation) -> npt.NDArray[np.float64]:
    "get separation/radial bin averages from allcounts or xi estimator"
    return allcount_switch_function(allcounts, lambda x: x.sepavg(axis=0), lambda x: x.get('RR').values('s'))


def get_mu_avg_from_allcounts(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation) -> npt.NDArray[np.float64]:
    "get mu/angular bin averages from allcounts or xi estimator"
    return allcount_switch_function(allcounts, lambda x: x.sepavg(axis=1), lambda x: x.get('RR').values('mu'))


def fix_and_wrap_pycorr(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator) -> pycorr.twopoint_estimator.BaseTwoPointEstimator:
    "check for negative counts in pycorr estimator, try to fix by symmetry and wrap to |mu| bins"
    return fix_bad_bins_pycorr(allcounts).wrap()


def fix_and_wrap_lsstypes(allcounts: lsstypes.Count2Correlation) -> lsstypes.Count2Correlation:
    "check for negative counts in lsstypes Count2Correlation, try to fix by symmetry and wrap to |mu| bins"
    return wrap_correlation(fix_bad_bins_lsstypes(allcounts))


def fix_and_wrap_allcounts(allcounts: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation) -> pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation:
    """
    For estimators with negative mu values, check for negative counts and wrap to |mu| bins for either pycorr estimator or lsstypes Count2Correlation.
    For estimators with only positive mu values, return the input as is, assuming it is already in |mu| bins.
    """
    return allcount_switch_function(allcounts, fix_and_wrap_pycorr, fix_and_wrap_lsstypes) if np.any(get_mu_edges_from_allcounts(allcounts) < 0) else allcounts # only try to fix and wrap if the first mu edge is negative, as a heuristic for whether it is already in |mu| bins