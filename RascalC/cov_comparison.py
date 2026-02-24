"Compact measures for covariance matrix comparison introduced in Section 3.1 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ and briefly explained again in Section 2.2 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_."
import numpy as np


def rms_eig_inv_test_covs(C1: np.typing.NDArray[np.float64], C2: np.typing.NDArray[np.float64]) -> float:
    "Compute the R_inv comparison measure (see Equation 3.4 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ or Equation 2.17 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_) between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    N = len(C2)
    tmp = Psi1.dot(C2) - np.eye(N)
    return np.sqrt(np.sum(np.diag(tmp.dot(tmp))) / N)


def KL_div_covs(C1: np.typing.NDArray[np.float64], C2: np.typing.NDArray[np.float64]) -> float:
    "Compute the Kullback-Leibler divergence between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return (np.trace(Psi1C2) - len(C2) - np.log(np.linalg.det(Psi1C2)))/2


def chi2_red_covs(C1: np.typing.NDArray[np.float64], C2: np.typing.NDArray[np.float64]) -> float:
    "Compute the reduced chi-squared comparison measure (see Equation 3.4 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ or Equation 2.18 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_) between two covariance matrices; the first is inverted."
    Psi1 = np.linalg.inv(C1)
    Psi1C2 = Psi1.dot(C2)
    return np.trace(Psi1C2)/len(C2)