## Script to post-process the multi-field integrals computed by the C++ code. This computes two shot-noise rescaling parameters, alphas, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_legendre, cov_filter_legendre_pycorr, load_matrices_multi, check_eigval_convergence, fit_shot_noise_rescaling, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_legendre
from typing import Iterable, Callable


def post_process_legendre_mocks_multi(mock_cov_file: str, file_root: str, n: int, max_l: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print) -> dict[str]:
    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)

    mock_cov = np.loadtxt(mock_cov_file)[cov_filter_legendre_pycorr(n, max_l, skip_r_bins, skip_l, multi=True)] # load external mock covariance matrix; select bins and switch their ordering right away

    input_file = load_raw_covariances_legendre(file_root, n, max_l, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load autocovariance from mock covariance
    n_bins = len(mock_cov) // 3
    mock_cov_11 = mock_cov[:n_bins, :n_bins]
    mock_cov_22 = mock_cov[2*n_bins:, 2*n_bins:]
    auto_mock_cov = [mock_cov_11, mock_cov_22]

    alpha_best = np.ones(2) # fill with ones by default, although this should not matter

    # Load full matrices
    c2, c3, c4 = load_matrices_multi(input_file, cov_filter, full = True, jack = False)
    # Load subsample matrices
    c2s, c3s, c4s = load_matrices_multi(input_file, cov_filter, full = False, jack = False)

    ## Optimize for alpha_1 and alpha_2 separately using single tracer auto-covariances
    for t, this_mock_cov in enumerate(auto_mock_cov):

        # Load single-tracer parts
        this_c2 = c2[t, t]
        this_c2s = c2s[t, t]
        this_c3 = c3[t, t, t]
        this_c3s = c3s[t, t, t]
        this_c4 = c4[t, t, t, t]
        this_c4s = c4s[t, t, t, t]

        # Check matrix convergence
        eigval_ok = check_eigval_convergence(this_c2, this_c4, f"Tracer {t+1}", print_function = print_function)

        # Now optimize for shot-noise rescaling parameter alpha
        print_function("Optimizing for the shot-noise rescaling parameter alpha_%d" % (t+1))
        optimal_alpha = fit_shot_noise_rescaling(this_mock_cov, this_c2, this_c3, this_c4, this_c2s, this_c3s, this_c4s)
        print_function("Optimization complete for field %d - optimal rescaling parameter is alpha_%d = %.6f" % (t+1, t+1, optimal_alpha))

        alpha_best[t] = optimal_alpha

        # Check matrix convergence for the optimal alpha: if it is <1, the eigenvalue criterion should be strengthened
        if eigval_ok and optimal_alpha < 1: check_eigval_convergence(this_c2, this_c4, optimal_alpha, kind = f"Tracer {t+1}")

    # Compute full matrices
    c_tot, c_comb = add_cov_terms_multi(c2, c3, c4, alpha_best)

    # Check positive definiteness
    check_positive_definiteness(c_comb)

    # Compute subsampled matrices (all submatrices combined)
    _, c_comb_subsamples = add_cov_terms_multi(c2s, c3s, c4s, alpha_best)

    # Now compute precision matrix
    D_est, prec_comb = compute_D_precision_matrix(c_comb_subsamples, c_comb)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff = compute_N_eff_D(D_est, print_function)

    output_dict = {"full_theory_covariance": c_comb, "all_covariances": c_tot, "shot_noise_rescaling": alpha_best, "full_theory_precision": prec_comb, "N_eff": N_eff, "full_theory_D_matrix": D_est, "individual_theory_covariances": c_comb_subsamples, "mock_covariance": mock_cov}

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Legendre_Mocks_n%d_m%d.npz' % (n, max_l))
    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s" % output_name)

    return output_dict
