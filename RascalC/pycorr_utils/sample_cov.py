r"These functions generate sample covariances of binned :math:`\xi(s,\mu)` from ``cosmodesi/pycorr`` ``s_mu`` `2PCF estimators <https://github.com/cosmodesi/pycorr>`_"

import pycorr
import numpy as np
from .utils import reshape_pycorr


def sample_cov_from_pycorr(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf) -> np.typing.NDArray[np.float64]:
    r"""
    Produce a sample covariance of binned :math:`\xi(s,\mu)` from ``cosmodesi/pycorr`` ``s_mu`` `2PCF estimators <https://github.com/cosmodesi/pycorr>`_.
    Multiple tracers are supported.

    Parameters
    ----------
    xi_estimators : list of lists of ``pycorr.twopoint_estimator.BaseTwoPointEstimator``\s
        The first element must be the list of first tracer auto-correlation function estimators.
        The (optional) second element should be the list of cross-correlation function estimators between the 1st and the 2nd tracer. (Ordering of mock realizations must be the same, or the cross-covariance blocks will be wrong.)
        The (optional) third element should be the list of second tracer auto-correlation function estimators.
        This ordering can be extended (or altered) but it must be consistent.
    
    n_mu : integer or None
        (Optional) If integer, number of angular (:math:`\mu`) bins to rebin to (after wrapping).
        If None (default), the angular bins are unchanged aside from wrapping of :math:`-1 \le \mu < 0` bins into :math:`0 \le \mu \le 1`.
    
    r_step: float or None
        (Optional) If float, sets the uniform spacing of radial/separation bins.
        If None (default), the radial/separation binning is not changed.
    
    r_max: float
        (Optional) Sets the maximum radius/separation, any bins beyond that value are discarded. By default, the value is infinity, so no bins are discarded.
    """
    if len(xi_estimators) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(xi_estimators[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(xi_estimators_c) != len(xi_estimators[0]) for xi_estimators_c in xi_estimators[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    # convert each xi estimator to binned xi array, and then turn list of lists into array too
    xi = np.array([[reshape_pycorr(xi_estimator, n_mu = n_mu, r_step = r_step, r_max = r_max).corr for xi_estimator in xi_estimators_c] for xi_estimators_c in xi_estimators])
    # now indices are [c, s, r, m]: correlation number, sample number, radial bin and then mu bin
    # need [s, c, r, m]
    xi = xi.transpose(1, 0, 2, 3)
    # now flatten all dimensions except the samples
    xi = xi.reshape(xi.shape[0], -1)
    return np.cov(xi.T) # xi has to be transposed, because variables (bins) are in columns (2nd index) of it and np.cov expects otherwise.
    # Weights are assumed the same, hard to figure out alternatives, and they do not seem necessary


def sample_cov_from_pycorr_to_file(xi_estimators: list[list[pycorr.twopoint_estimator.BaseTwoPointEstimator]], outfile_name: str, n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf) -> None:
    r"""
    Produce a sample covariance of binned :math:`\xi(s,\mu)` from ``cosmodesi/pycorr`` ``s_mu`` `2PCF estimators <https://github.com/cosmodesi/pycorr>`_ and write the matrix to a text file.
    Multiple tracers are supported.

    Parameters are analogous to :func:`sample_cov_from_pycorr` except
    -----------------------------------------------------------------
    outfile_name : string (filename)
        The name for the output text file.
    """
    np.savetxt(outfile_name, sample_cov_from_pycorr(xi_estimators, n_mu, r_step, r_max))


def sample_cov_from_pycorr_files(infile_names: list[list[str]], outfile_name: str, n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf):
    r"""
    Produce a sample covariance of binned :math:`\xi(s,\mu)` from ``cosmodesi/pycorr`` ``s_mu`` `.npy` files and write the matrix to a text file.
    Multiple tracers are supported.

    Parameters
    ----------
    infile_names : list of lists of strings (filenames)
        The first element must be the list of first tracer auto-correlation function estimator filenames.
        The (optional) second element should be the list of cross-correlation function estimator filenames between the 1st and the 2nd tracer. (Ordering of mock realizations must be the same, or the cross-covariance blocks will be wrong.)
        The (optional) third element should be the list of second tracer auto-correlation function estimator filenames.
        This ordering can be extended (or altered) but it must be consistent.

    outfile_name : string
        The name for the output text file.
    
    n_mu : integer or None
        (Optional) If integer, number of angular (:math:`\mu`) bins to rebin to (after wrapping).
        If None (default), the angular bins are unchanged aside from wrapping of :math:`-1 \le \mu < 0` bins into :math:`0 \le \mu \le 1`.
    
    r_step: float or None
        (Optional) If float, sets the uniform spacing of radial/separation bins.
        If None (default), the radial/separation binning is not changed.
    
    r_max: float
        (Optional) Sets the maximum radius/separation, any bins beyond that value are discarded. By default, the value is infinity, so no bins are discarded.
    """
    if len(infile_names) <= 0: raise ValueError("Need at least one correlation function group in the outer list")
    if len(infile_names[0]) < 2: raise ValueError("Need at least two samples to compute the covariance matrix")
    if any(len(infile_names_c) != len(infile_names[0]) for infile_names_c in infile_names[1:]):
        raise ValueError("Need the same number of files for different correlation functions")
    xi_estimators = [[reshape_pycorr(pycorr.TwoPointCorrelationFunction.load(infile_name), n_mu = n_mu, r_step = r_step, r_max = r_max) for infile_name in infile_names_c] for infile_names_c in infile_names]
    sample_cov_from_pycorr_to_file(xi_estimators, outfile_name, n_mu, r_step, r_max)
