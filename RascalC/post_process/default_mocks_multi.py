## Script to post-process the multi-field integrals computed by the C++ code. This computes two shot-noise rescaling parameters, alphas, from a mock derived covariance matrix.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_smu, load_matrices_multi, check_eigval_convergence, fit_shot_noise_rescaling, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_smu


def post_process_default_mocks_multi(mock_cov_file: str, file_root: str, n: int, m: int, outdir: str, skip_r_bins: int = 0, print_function = print) -> dict[str]:
    skip_bins = skip_r_bins * m
    n_bins = n * m - skip_bins

    skip_mask = np.tile(np.arange(n * m) >= skip_bins, 3) # the mask gives False for first skip_bins, all this repeating 3 times; 3 is number of correlations for 2 tracers
    mock_cov = np.loadtxt(mock_cov_file)[np.ix_(skip_mask, skip_mask)] # load external mock covariance matrix, select bins on both axes

    cov_filter = cov_filter_smu(n, m, skip_r_bins)

    input_file = load_raw_covariances_smu(file_root, n, m, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load autocovariance from mock covariance
    mock_cov_11 = mock_cov[:n_bins, :n_bins]
    mock_cov_22 = mock_cov[2*n_bins:, 2*n_bins:]
    auto_mock_cov = [mock_cov_11, mock_cov_22]

    alpha_best = np.zeros(2)

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
        check_eigval_convergence(this_c2, this_c4, f"Tracer {t+1}")

        # Now optimize for shot-noise rescaling parameter alpha
        print_function("Optimizing for the shot-noise rescaling parameter alpha_%d" % (t+1))
        optimal_alpha = fit_shot_noise_rescaling(this_mock_cov, this_c2, this_c3, this_c4, this_c2s, this_c3s, this_c4s)
        print_function("Optimization complete for field %d - optimal rescaling parameter is alpha_%d = %.6f" % (t+1, t+1, optimal_alpha))

        alpha_best[t] = optimal_alpha

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

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Default_Mocks_n%d_m%d.npz' % (n, m))
    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s" % output_name)

    return output_dict
