## Script to post-process the single-field integrals computed by the C++ code in mixed Legendre (LEGENDRE_MIX) mode. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import numpy as np
import sys, os
from warnings import warn
from .utils import cov_filter_legendre, load_matrices_single, check_eigval_convergence, add_cov_terms_single, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D, fit_shot_noise_rescaling
from .collect_raw_covariance_matrices import load_raw_covariances_legendre


def post_process_legendre_mix_jackknife(jackknife_file: str, weight_dir: str, file_root: str, m: int, max_l: int, outdir: str, skip_r_bins: int = 0, skip_l: int = 0, tracer: int = 1, print_function = print) -> dict[str]:
    # Load jackknife xi estimates from data
    print_function("Loading correlation function jackknife estimates from %s" % jackknife_file)
    xi_jack = np.loadtxt(jackknife_file, skiprows = 2)
    n_jack = xi_jack.shape[0] # total jackknives
    n = xi_jack.shape[1] // m # radial bins
    n_bins = (n_l - skip_l) * (n - skip_r_bins) # total Legendre bins to work with

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
    data_cov = data_cov[skip_r_bins:, :, skip_r_bins:, :] # discard the extra radial bins now since it is convenient
    data_cov = np.einsum("imjn,mp,nq->ipjq", data_cov, mu_bin_legendre_factors, mu_bin_legendre_factors) # use mu bin Legendre factors to project mu bins into Legendre multipoles, staying within the same radial bins. The indices are now [r_bin, ell] for rows and columns
    data_cov = data_cov.reshape(n_bins, n_bins)

    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)
    n_l = max_l // 2 + 1 # number of multipoles
    
    input_file = load_raw_covariances_legendre(file_root, n, max_l, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full jackknife theoretical matrices
    print_function("Loading best estimate of jackknife covariance matrix")
    c2j, c3j, c4j = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = True)

    # Check matrix convergence
    check_eigval_convergence(c2j, c4j, "Jackknife")

    # Load in partial jackknife theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = True)

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    alpha_best = fit_shot_noise_rescaling(data_cov, c2j, c3j, c4j, c2s, c3s, c4s)
    print_function("Optimization complete - optimal rescaling parameter is %.6f" % alpha_best)

    # Compute jackknife and full covariance matrices
    jack_cov = add_cov_terms_single(c2j, c3j, c4j, alpha_best)
    partial_jack_cov = add_cov_terms_single(c2s, c3s, c4s, alpha_best)
    _, jack_prec = compute_D_precision_matrix(partial_jack_cov, jack_cov)

    c2f, c3f, c4f = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)
    full_cov = add_cov_terms_single(c2f, c3f, c4f, alpha_best)

    # Check convergence
    check_eigval_convergence(c2f, c4f, "Full")

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

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_Jackknife_n%d_l%d_j%d.npz' % (n, max_l, n_jack))
    np.savez(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (7, 8, 9):
        print("Usage: python post_process_legendre_mix_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {MAX_L} {OUTPUT_DIR} [{SKIP_R_BINS} [{SKIP_L}]]")
        sys.exit(1)
            
    jackknife_file = str(sys.argv[1])
    weight_dir = str(sys.argv[2])
    file_root = str(sys.argv[3])
    m = int(sys.argv[4])
    max_l = int(sys.argv[5])
    outdir = str(sys.argv[6])
    from .utils import get_arg_safe
    skip_r_bins = get_arg_safe(7, int, 0)
    skip_l = get_arg_safe(8, int, 0)

    post_process_legendre_mix_jackknife(jackknife_file, weight_dir, file_root, m, max_l, outdir, skip_r_bins, skip_l)