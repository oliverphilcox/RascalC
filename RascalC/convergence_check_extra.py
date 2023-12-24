## Script to perform an extra convergence check on full integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Should work in any case - default, jackknife, Legendre, multi-tracer - as it utilizes universal data from RascalC file

import numpy as np
from .utils import blank_function


# methods to assess similarity
def rms_eig_inv_test_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    Psi1 = np.linalg.inv(C1)
    N = len(C2)
    tmp = Psi1.dot(C2) - np.eye(N)
    return np.sqrt(np.sum(np.diag(tmp.dot(tmp))) / N)


def KL_div_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return (np.trace(Psi1C2) - len(C2) - np.log(np.linalg.det(Psi1C2)))/2


def chi2_red_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return np.trace(Psi1C2)/len(C2)


def cmp_cov(cov_first: np.ndarray[float], cov_second: np.ndarray[float], print_function = blank_function) -> dict[str, float]:
    result = dict()

    result["R_inv"] = (rms_eig_inv_test_covs(cov_first, cov_second), rms_eig_inv_test_covs(cov_second, cov_first))
    print_function("RMS eigenvalues of inverse tests for cov half-estimates are %.2e and %.2e" % result["R_inv"])

    result["D_KL"] = (KL_div_covs(cov_first, cov_second), KL_div_covs(cov_second, cov_first))
    print_function("KL divergences between cov half-estimates are %.2e and %.2e" % result["D_KL"])

    result["chi2_red-1"] = (chi2_red_covs(cov_first, cov_second)-1, chi2_red_covs(cov_second, cov_first)-1)
    print_function("Reduced chi2-1 between cov half-estimates are %.2e and %.2e" % result["chi2_red-1"])

    return result


def convergence_check_extra_splittings(c_samples: np.ndarray[float], n_samples: int | None = None, print_function = blank_function) -> dict[str, dict[str, float]]:
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


def convergence_check_extra(rascalc_results: dict[str], n_samples: int | None = None, print_function = blank_function) -> dict[str, dict[str, dict[str, float]]]:
    print_function("Full covariance")
    result = {"full": convergence_check_extra_splittings(rascalc_results["individual_theory_covariances"], n_samples, print_function)}

    jack_key = "individual_theory_jackknife_covariances"
    if jack_key in rascalc_results.keys():
        print_function("Jack covariance")
        result["jack"] = convergence_check_extra_splittings(rascalc_results[jack_key], n_samples, print_function)
    return result


def convergence_check_extra_file(rascalc_results_file: str, n_samples: int | None = None, print_function = blank_function) -> dict[str, dict[str, dict[str, float]]]:
    with np.load(rascalc_results_file) as f:
        return convergence_check_extra(f, n_samples, print_function)