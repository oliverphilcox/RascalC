## Script to post-process the single-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_smu, load_matrices_single, check_eigval_convergence, add_cov_terms_single, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D, fit_shot_noise_rescaling
from ..raw_covariance_matrices import load_raw_covariances_smu
from typing import Literal, Callable, Iterable


def post_process_default_mocks(mock_cov_file: str, file_root: str, n: int, m: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, tracer: Literal[1, 2] = 1, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print) -> dict[str]:
    cov_filter = cov_filter_smu(n, m, skip_r_bins)
    mock_cov = np.loadtxt(mock_cov_file)[cov_filter] # load external mock covariance matrix

    input_file = load_raw_covariances_smu(file_root, n, m, n_samples, print_function)

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

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Default_Mocks_n%d_m%d.npz'%(n,m))
    np.savez_compressed(output_name, **output_dict)
    output_dict["path"] = output_name
    output_dict["filename"] = os.path.basename(output_name)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict
