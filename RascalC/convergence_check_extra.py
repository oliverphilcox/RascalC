"""
Functions to perform extra convergence check on full (and jackknife) RascalC integrals.
More specifically, divide integral subsamples into halves and check similarity of their average results.
These methods work in any mode — e.g. jackknife, Legendre, multi-tracer — as it utilizes universal data from RascalC file
"""

import numpy as np
from .utils import blank_function
from .cov_comparison import rms_eig_inv_test_covs, KL_div_covs, chi2_red_covs
from typing import Callable


def cmp_cov(cov_first: np.ndarray[float], cov_second: np.ndarray[float], print_function: Callable[[str], None] = blank_function) -> dict[str, float]:
    """
    Compute the selected comparison measures between two covariance matrices and return as a dictionary.
    This method is decribed in Section 3.2 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_.
    Optionally, use ``print_function`` to report the results.
    """
    result = dict()

    result["R_inv"] = (rms_eig_inv_test_covs(cov_first, cov_second), rms_eig_inv_test_covs(cov_second, cov_first))
    print_function("R_inv (RMS eigenvalues of inverse tests) for cov half-estimates are %.2e and %.2e" % result["R_inv"])

    result["D_KL"] = (KL_div_covs(cov_first, cov_second), KL_div_covs(cov_second, cov_first))
    print_function("KL divergences between cov half-estimates are %.2e and %.2e" % result["D_KL"])

    result["chi2_red-1"] = (chi2_red_covs(cov_first, cov_second)-1, chi2_red_covs(cov_second, cov_first)-1)
    print_function("Reduced chi2-1 between cov half-estimates are %.2e and %.2e" % result["chi2_red-1"])

    return result


def convergence_check_extra_splittings(c_samples: np.ndarray[float], n_samples: int | None = None, print_function: Callable[[str], None] = blank_function) -> dict[str, dict[str, float]]:
    """
    Perform two different splittings in halves using the covariance matrix samples ``c_samples``, compute the comparison measures between the two average covariance matrices and return as a dictionary.
    This method is decribed in Section 3.2 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_.
    Optionally, use only ``n_samples`` first samples.
    Further optionally, use ``print_function`` to report the results.
    """
    if n_samples is None: n_samples = len(c_samples)
    n_samples_2 = n_samples // 2

    result = dict()

    print_function("First splitting")
    cov_first = np.mean(c_samples[:n_samples_2], axis=0)
    cov_second = np.mean(c_samples[n_samples_2:n_samples], axis=0)
    result["split1"] = cmp_cov(cov_first, cov_second, print_function)

    print_function("Second splitting")
    cov_first = np.mean(c_samples[:n_samples:2], axis=0)
    cov_second = np.mean(c_samples[1:n_samples:2], axis=0)
    result["split2"] = cmp_cov(cov_first, cov_second, print_function)

    return result


def convergence_check_extra(rascalc_results: dict[str], n_samples: int | None = None, print_function: Callable[[str], None] = blank_function) -> dict[str, dict[str, dict[str, float]]]:
    """
    Perform two different splittings in halves using the RascalC results file/dictionary, compute the comparison measures between the two average covariance matrices and return as a dictionary.
    Do this for full and jackknife covariance matrices (if the latter are present).
    This method is decribed in Section 3.2 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_.
    Optionally, use only ``n_samples`` first samples.
    Further optionally, use ``print_function`` to report the results.
    """
    print_function("Full covariance")
    result = {"full": convergence_check_extra_splittings(rascalc_results["individual_theory_covariances"], n_samples, print_function)}

    jack_key = "individual_theory_jackknife_covariances"
    if jack_key in rascalc_results.keys():
        print_function("Jack covariance")
        result["jack"] = convergence_check_extra_splittings(rascalc_results[jack_key], n_samples, print_function)
    return result


def convergence_check_extra_file(rascalc_results_filename: str, n_samples: int | None = None, print_function: Callable[[str], None] = blank_function) -> dict[str, dict[str, dict[str, float]]]:
    """
    Perform two different splittings in halves using the RascalC results filename, compute the comparison measures between the two average covariance matrices and return as a dictionary.
    Do this for full and jackknife covariance matrices (if the latter are present).
    This method is decribed in Section 3.2 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_.
    Optionally, use only ``n_samples`` first samples.
    Further optionally, use ``print_function`` to report the results.
    """
    with np.load(rascalc_results_filename) as f:
        return convergence_check_extra(f, n_samples, print_function)