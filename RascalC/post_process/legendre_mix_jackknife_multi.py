## Function to post-process the multi-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data-derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import numpy as np
import os
from warnings import warn
from .utils import format_skip_r_bins, cov_filter_legendre, load_matrices_multi, check_eigval_convergence, fit_shot_noise_rescaling, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_legendre
from typing import Iterable, Callable


def post_process_legendre_mix_jackknife_multi(jackknife_file_11: str, jackknife_file_12: str, jackknife_file_22: str, weight_dir: str, file_root: str, m: int, max_l: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print):
    ## First load jackknife xi estimates from data:
    print_function("Loading correlation function jackknife estimates")
    xi_jack_11 = np.loadtxt(jackknife_file_11, skiprows=2)
    xi_jack_12 = np.loadtxt(jackknife_file_12, skiprows=2)
    xi_jack_22 = np.loadtxt(jackknife_file_22, skiprows=2)
    if not (xi_jack_11.shape == xi_jack_22.shape == xi_jack_12.shape): raise ValueError('Must have the same configuration of jackknives for each field.')

    n_bins_smu = xi_jack_11.shape[1] # total bins
    n_jack = xi_jack_11.shape[0] # total jackknives
    n = n_bins_smu // m # radial bins
    n_l = max_l // 2 + 1 # number of even multipoles
    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    n_bins = (n_l - skip_l) * (n - skip_r_bins_start - skip_r_bins_end) # total Legendre bins to work with

    # First exclude any dodgy jackknife regions
    good_jk = np.where(np.all(np.isfinite(xi_jack_11) & np.isfinite(xi_jack_12) & np.isfinite(xi_jack_22), axis=1))[0] # all xi in jackknife have to be normal numbers
    if len(good_jk) < n_jack:
        warn("Using only %d out of %d jackknives - some xi values were not finite" % (len(good_jk), n_jack))

    xi_all = np.zeros([len(good_jk), 3*n_bins_smu])
    weights_all = np.zeros_like(xi_all)

    # Load in all xi functions
    xi_all[:, :n_bins_smu] = xi_jack_11[good_jk]
    xi_all[:, n_bins_smu:2*n_bins_smu] = xi_jack_12[good_jk]
    xi_all[:, 2*n_bins_smu:] = xi_jack_22[good_jk]

    # Load in all weights:
    weight_file11 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack))
    weight_file12 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_12.dat'%(n,m,n_jack))
    weight_file22 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_22.dat'%(n,m,n_jack))
    weights11 = np.loadtxt(weight_file11)[:, 1:]
    weights12 = np.loadtxt(weight_file12)[:, 1:]
    weights22 = np.loadtxt(weight_file22)[:, 1:]
    weights = np.array([these_weights[good_jk] for these_weights in (weights11, weights12, weights22)])
    weights /= np.sum(weights, axis = 1) # renormalize after possibly discarding some jackknives
    weights_all = weights.swapaxes(0, 1).reshape(len(good_jk), 3*n_bins_smu)

    # Compute full covariance matrix:
    tmp_cov = np.zeros([len(good_jk), 3*n_bins_smu])
    mean_xi = np.sum(xi_all * weights_all, axis=0)
    tmp_cov = weights_all * (xi_all-mean_xi)

    print_function("Computing full data covariance matrix")
    # Now compute covariance matrix:
    num = np.matmul(tmp_cov.T, tmp_cov)
    denom = np.matmul(weights_all.T, weights_all)
    data_cov = num/(np.ones_like(denom)-denom)

    mu_bin_legendre_file = os.path.join(weight_dir, 'mu_bin_legendre_factors_m%d_l%d.txt' % (m, max_l))
    print("Loading mu bin Legendre factors from %s" % mu_bin_legendre_file)
    mu_bin_legendre_factors = np.loadtxt(mu_bin_legendre_file) # rows correspond to mu bins, columns to multipoles
    if skip_l > 0: mu_bin_legendre_factors = mu_bin_legendre_factors[:, :-skip_l] # discard unneeded l; the expression works wrong for skip_l=0

    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins)

    input_file = load_raw_covariances_legendre(file_root, n, max_l, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load autocovariance from data-covariance
    data_cov_11 = data_cov[:n_bins_smu,:n_bins_smu]
    data_cov_22 = data_cov[2*n_bins_smu:,2*n_bins_smu:]
    auto_data_cov = [data_cov_11,data_cov_22]

    alpha_best = np.ones(2) # fill with ones by default, although this should not matter

    # Load full jack matrices
    c2j, c3j, c4j = load_matrices_multi(input_file, cov_filter, full = True, jack = True)
    # Load subsample jack matrices
    c2s, c3s, c4s = load_matrices_multi(input_file, cov_filter, full = False, jack = True)

    ## Optimize for alpha_1 and alpha_2 separately.
    for t, this_data_cov in enumerate(auto_data_cov):

        # Project the data jackknife covariance from mu bins to Legendre multipoles
        this_data_cov = this_data_cov.reshape(n, m, n, m) # make the array 4D with [r_bin, mu_bin] indices for rows and columns
        this_data_cov = this_data_cov[skip_r_bins_start:, :, skip_r_bins_start:, :] # discard the extra radial bins now since it is convenient
        if skip_r_bins_end != 0: this_data_cov = this_data_cov[:-skip_r_bins_end, :, :-skip_r_bins_end, :] # discard radial bins at the end if any
        this_data_cov = np.einsum("imjn,mp,nq->ipjq", this_data_cov, mu_bin_legendre_factors, mu_bin_legendre_factors) # use mu bin Legendre factors to project mu bins into Legendre multipoles, staying within the same radial bins. The indices are now [r_bin, ell] for rows and columns
        this_data_cov = this_data_cov.reshape(n_bins, n_bins)

        # Load single-tracer parts
        this_c2j = c2j[t, t]
        this_c2s = c2s[t, t]
        this_c3j = c3j[t, t, t]
        this_c3s = c3s[t, t, t]
        this_c4j = c4j[t, t, t, t]
        this_c4s = c4s[t, t, t, t]

        # Check matrix convergence
        eigval_ok = check_eigval_convergence(this_c2j, this_c4j, kind = f"Tracer {t+1} jackknife", print_function = print_function)

        # Now optimize for shot-noise rescaling parameter alpha
        print_function("Optimizing for the shot-noise rescaling parameter alpha_%d" % (t+1))
        optimal_alpha = fit_shot_noise_rescaling(this_data_cov, this_c2j, this_c3j, this_c4j, this_c2s, this_c3s, this_c4s)
        print_function("Optimization complete for field %d - optimal rescaling parameter is alpha_%d = %.6f" % (t+1, t+1,optimal_alpha))

        alpha_best[t] = optimal_alpha

        # Check matrix convergence for the optimal alpha: if it is <1, the eigenvalue criterion should be strengthened
        if eigval_ok and optimal_alpha < 1: check_eigval_convergence(this_c2j, this_c4j, optimal_alpha, kind = f"Tracer {t+1} jackknife")

    # Load full matrices
    c2f, c3f, c4f = load_matrices_multi(input_file, cov_filter, full = True, jack = False)
    c_tot, c_comb = add_cov_terms_multi(c2f, c3f, c4f, alpha_best)

    # Check positive definiteness
    check_positive_definiteness(c_comb)

    # Compute subsampled matrices (all submatrices combined)
    c2fs, c3fs, c4fs = load_matrices_multi(input_file, cov_filter, full = False, jack = False)
    _, c_comb_subsamples = add_cov_terms_multi(c2fs, c3fs, c4fs, alpha_best)
    
    # Compute jackknive totals
    cj_tot, cj_comb = add_cov_terms_multi(c2j, c3j, c4j, alpha_best)
    _, cj_comb_subsamples = add_cov_terms_multi(c2s, c3s, c4s, alpha_best)

    # Now compute precision matrices
    D_est, prec_comb = compute_D_precision_matrix(c_comb_subsamples, c_comb)
    Dj_est, precj_comb = compute_D_precision_matrix(cj_comb_subsamples, cj_comb)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff = compute_N_eff_D(D_est, print_function)
    # Nj_eff = compute_N_eff_D(Dj_est, print_function)

    output_dict = {"full_theory_covariance": c_comb, "all_covariances": c_tot, "shot_noise_rescaling": alpha_best, "full_theory_precision": prec_comb, "N_eff": N_eff, "full_theory_D_matrix": D_est, "individual_theory_covariances": c_comb_subsamples, "jackknife_data_covariance": data_cov, "jackknife_theory_covariance": cj_comb, "all_jackknife_covariances": cj_tot, "jackknife_theory_precision": precj_comb, "individual_theory_jackknife_covariances": cj_comb_subsamples}

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Legendre_Jackknife_n%d_m%d_j%d.npz' % (n, max_l, n_jack))

    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s" % output_name)

    return output_dict
