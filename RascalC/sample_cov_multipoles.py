r"These functions generate sample covariances of binned :math:`\xi_\ell(s)` from ``s_mu`` 2PCF estimators from `cosmodesi/pycorr <https://github.com/cosmodesi/pycorr>`_ or `adematti/lsstypes <https://github.com/adematti/lsstypes>`_"

import pycorr
import lsstypes
import numpy as np
import numpy.typing as npt
from .allcounts_utils import allcount_switch_function


def get_poles(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation, ells: list[int]) -> np.ndarray:
    "get binned correlation function multipolesarray from either pycorr estimator or lsstypes Count2Correlation"
    return allcount_switch_function(xi_estimator, lambda x: x.get_corr(mode='poles', ells=ells), lambda x: x.project(mode='poles', ells=ells))


def sample_cov_multipoles(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation]], max_l: int) -> npt.NDArray[np.float64]:
    r"""
    Produce a sample covariance of binned :math:`\xi_\ell(s)` from ``s_mu`` 2PCF estimators from `cosmodesi/pycorr <https://github.com/cosmodesi/pycorr>`_ or `adematti/lsstypes <https://github.com/adematti/lsstypes>`_.
    Considers only even multipoles.
    Multiple tracers are supported.

    Parameters
    ----------
    xi_estimators : list of lists of :class:`pycorr.twopoint_estimator.BaseTwoPointEstimator`\s or :class:`lsstypes.Count2Correlation`\s
        The first element must be the list of first tracer auto-correlation function estimators.
        The (optional) second element should be the list of cross-correlation function estimators between the 1st and the 2nd tracer. (Ordering of mock realizations must be the same, or the cross-covariance blocks will be wrong.)
        The (optional) third element should be the list of second tracer auto-correlation function estimators.
        This ordering can be extended (or altered) but it must be consistent.
    
    max_l : integer
        Max (even) Legendre multipole index.
    """
    if max_l < 0: raise ValueError("Maximal multipole can not be negative")
    assert max_l % 2 == 0, "Only even multipoles supported"
    if len(xi_estimators) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(xi_estimators[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(xi_estimators_c) != len(xi_estimators[0]) for xi_estimators_c in xi_estimators[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    ells = np.arange(0, max_l+1, 2)
    # convert each xi estimator to multipoles array, and then turn list of lists into array too
    xi = np.array([[get_poles(xi_estimator, ells) for xi_estimator in xi_estimators_c] for xi_estimators_c in xi_estimators])
    # now indices are [c, s, l, r]: correlation number, sample number, multipole index and then radial bin
    # need [s, c, l, r]
    xi = xi.transpose(1, 0, 2, 3)
    # now flatten all dimensions except the samples
    xi = xi.reshape(xi.shape[0], -1)
    return np.cov(xi.T) # xi has to be transposed, because variables (bins) are in columns (2nd index) of it and np.cov expects otherwise.
    # Weights are assumed the same, hard to figure out alternatives, and they do not seem necessary


def sample_cov_multipoles_to_file(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator | lsstypes.Count2Correlation]], outfile_name: str, max_l: int) -> None:
    r"""
    Produce a sample covariance of binned :math:`\xi_\ell(s)` from ``s_mu`` 2PCF estimators from `cosmodesi/pycorr <https://github.com/cosmodesi/pycorr>`_ or `adematti/lsstypes <https://github.com/adematti/lsstypes>`_ and write the matrix to a text file.
    Considers only even multipoles.
    Multiple tracers are supported.

    Parameters are analogous to :func:`sample_cov_multipoles` except
    ----------------------------------------------------------------------------
    outfile_name : string (filename)
        The name for the output text file.
    """
    np.savetxt(outfile_name, sample_cov_multipoles(xi_estimators, max_l))
