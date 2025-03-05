"Compact measures for covariance matrix comparison introduced in `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ and briefly explained again in `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_."
import numpy as np


def rms_eig_inv_test_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    "Compute the R_inv comparison measure (see `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ or `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_) between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    N = len(C2)
    tmp = Psi1.dot(C2) - np.eye(N)
    return np.sqrt(np.sum(np.diag(tmp.dot(tmp))) / N)


def KL_div_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    "Compute the Kullback-Leibler divergence between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return (np.trace(Psi1C2) - len(C2) - np.log(np.linalg.det(Psi1C2)))/2


def chi2_red_covs(C1: np.ndarray[float], C2: np.ndarray[float]) -> float:
    "Compute the reduced chi-squared comparison measure (see `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ or `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_) between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return np.trace(Psi1C2)/len(C2)