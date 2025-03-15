## Script to post-process the single-field Legendre binned integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_legendre, load_matrices_single, check_eigval_convergence, add_cov_terms_single, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D, fit_shot_noise_rescaling
from ..raw_covariance_matrices import load_raw_covariances_legendre
from typing import Literal, Callable, Iterable


def cov_filter_legendre_mocks(n: int, max_l: int, skip_r_bins: int = 0, skip_l: int = 0):
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    l_indices = np.arange(n_l - skip_l)
    r_indices = np.arange(skip_r_bins, n)
    indices_l_r = (n * l_indices)[None, :] + r_indices[:, None] # this will transpose from pycorr to RascalC convention
    indices_1d = indices_l_r.ravel()
    return np.ix_(indices_1d, indices_1d)


def post_process_legendre_mocks(mock_cov_file: str, file_root: str, n: int, max_l: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, tracer: Literal[1, 2] = 1, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print) -> dict[str]:
    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)
    cov_filter_mocks = cov_filter_legendre_mocks(n, max_l, skip_r_bins, skip_l)

    mock_cov = np.loadtxt(mock_cov_file)[cov_filter_mocks] # load external mock covariance matrix
    
    input_file = load_raw_covariances_legendre(file_root, n, max_l, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c2f, c3f, c4f = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)

    # Check matrix convergence
    eigval_ok = check_eigval_convergence(c2f, c4f, print_function = print_function)

    # Load in partial theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = False)

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    alpha_best = fit_shot_noise_rescaling(mock_cov, c2f, c3f, c4f, c2s, c3s, c4s)
    print_function("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)

    # Check matrix convergence for the optimal alpha: if it is <1, the eigenvalue criterion should be strengthened
    if eigval_ok and alpha_best < 1: check_eigval_convergence(c2f, c4f, alpha_best)

    # Compute full covariance matrices and precision
    full_cov = add_cov_terms_single(c2f, c3f, c4f, alpha_best)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    partial_cov = add_cov_terms_single(c2s, c3s, c4s, alpha_best)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)

    output_dict = {"full_theory_covariance": full_cov, "shot_noise_rescaling": alpha_best, "full_theory_precision": full_prec, "N_eff": N_eff_D, "full_theory_D_matrix": full_D_est, "individual_theory_covariances": partial_cov, "mock_covariance": mock_cov}

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_Mocks_n%d_l%d.npz' % (n, max_l))
    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s" % output_name)

    return output_dict
