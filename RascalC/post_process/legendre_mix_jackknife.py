"""
Function to post-process the single-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.
"""

import numpy as np
import os
from warnings import warn
from .utils import format_skip_r_bins, cov_filter_legendre, load_matrices_single, check_eigval_convergence, add_cov_terms_single, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D, fit_shot_noise_rescaling
from ..raw_covariance_matrices import load_raw_covariances_legendre
from typing import Literal, Callable, Iterable


def post_process_legendre_mix_jackknife(jackknife_file: str, weight_dir: str, file_root: str, m: int, max_l: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, tracer: Literal[1, 2] = 1, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print, dry_run: bool = False) -> dict[str]:
    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_Jackknife_n%d_l%d_j%d.npz' % (n, max_l, n_jack))
    name_dict = dict(path=output_name, filename=os.path.basename(output_name))
    if dry_run: return name_dict

    # Load jackknife xi estimates from data
    print_function("Loading correlation function jackknife estimates from %s" % jackknife_file)
    xi_jack = np.loadtxt(jackknife_file, skiprows = 2)
    n_jack = xi_jack.shape[0] # total jackknives
    n = xi_jack.shape[1] // m # radial bins
    n_l = max_l // 2 + 1 # number of even multipoles
    skip_r_bins_start, skip_r_bins_end = format_skip_r_bins(skip_r_bins)
    n_bins = (n_l - skip_l) * (n - skip_r_bins_start - skip_r_bins_end) # total Legendre bins to work with

    weight_file = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_11.dat' % (n, m, n_jack))
    mu_bin_legendre_file = os.path.join(weight_dir, 'mu_bin_legendre_factors_m%d_l%d.txt' % (m, max_l))

    print_function("Loading jackknife weights from %s" % weight_file)
    weights = np.loadtxt(weight_file)[:, 1:]

    # First exclude any dodgy jackknife regions
    good_jk = np.where(np.all(np.isfinite(xi_jack), axis=1))[0] # all xi in jackknife have to be normal numbers
    if len(good_jk) < n_jack:
        warn("Using only %d out of %d jackknives - some xi values were not finite" % (len(good_jk), n_jack))
        xi_jack = xi_jack[good_jk]
        weights = weights[good_jk]
    weights /= np.sum(weights, axis=0) # renormalize weights after possibly discarding some jackknives

    # Compute data covariance matrix
    print_function("Computing data covariance matrix")
    mean_xi = np.sum(xi_jack * weights, axis = 0)
    tmp = weights * (xi_jack - mean_xi)
    data_cov = np.matmul(tmp.T, tmp)
    denom = np.matmul(weights.T, weights)
    data_cov /= (np.ones_like(denom) - denom)

    print("Loading mu bin Legendre factors from %s" % mu_bin_legendre_file)
    mu_bin_legendre_factors = np.loadtxt(mu_bin_legendre_file) # rows correspond to mu bins, columns to multipoles
    if skip_l > 0: mu_bin_legendre_factors = mu_bin_legendre_factors[:, :-skip_l] # discard unneeded l; the expression works wrong for skip_l=0

    # Project the data jackknife covariance from mu bins to Legendre multipoles
    data_cov = data_cov.reshape(n, m, n, m) # make the array 4D with [r_bin, mu_bin] indices for rows and columns
    data_cov = data_cov[skip_r_bins_start:, :, skip_r_bins_start:, :] # discard the extra radial bins now since it is convenient
    if skip_r_bins_end != 0: data_cov = data_cov[:-skip_r_bins_end, :, :-skip_r_bins_end, :] # discard radial bins at the end if any
    data_cov = np.einsum("imjn,mp,nq->ipjq", data_cov, mu_bin_legendre_factors, mu_bin_legendre_factors) # use mu bin Legendre factors to project mu bins into Legendre multipoles, staying within the same radial bins. The indices are now [r_bin, ell] for rows and columns
    data_cov = data_cov.reshape(n_bins, n_bins)

    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)
    n_l = max_l // 2 + 1 # number of multipoles
    
    input_file = load_raw_covariances_legendre(file_root, n, max_l, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full jackknife theoretical matrices
    print_function("Loading best estimate of jackknife covariance matrix")
    c2j, c3j, c4j = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = True)

    # Check matrix convergence
    eigval_ok = check_eigval_convergence(c2j, c4j, kind = "Jackknife", print_function = print_function)

    # Load in partial jackknife theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = True)

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    alpha_best = fit_shot_noise_rescaling(data_cov, c2j, c3j, c4j, c2s, c3s, c4s)
    print_function("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)

    # Check matrix convergence for the optimal alpha: if it is <1, the eigenvalue criterion should be strengthened
    if eigval_ok and alpha_best < 1: check_eigval_convergence(c2j, c4j, alpha_best, kind = "Jackknife")

    # Compute jackknife and full covariance matrices
    jack_cov = add_cov_terms_single(c2j, c3j, c4j, alpha_best)
    partial_jack_cov = add_cov_terms_single(c2s, c3s, c4s, alpha_best)
    _, jack_prec = compute_D_precision_matrix(partial_jack_cov, jack_cov)

    c2f, c3f, c4f = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)
    full_cov = add_cov_terms_single(c2f, c3f, c4f, alpha_best)

    # Check convergence
    check_eigval_convergence(c2f, c4f, alpha_best, kind = "Full", print_function = print_function)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial jackknife theoretical matrices
    c2fs, c3fs, c4fs = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = False)
    partial_cov = add_cov_terms_single(c2fs, c3fs, c4fs, alpha_best)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")    

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)

    output_dict = {"jackknife_theory_covariance": jack_cov, "full_theory_covariance": full_cov, "jackknife_data_covariance": data_cov, "shot_noise_rescaling": alpha_best, "jackknife_theory_precision": jack_prec, "full_theory_precision": full_prec, "N_eff": N_eff_D, "full_theory_D_matrix": full_D_est, "individual_theory_covariances": partial_cov, "individual_theory_jackknife_covariances": partial_jack_cov}

    np.savez_compressed(output_name, **output_dict)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict | name_dict
