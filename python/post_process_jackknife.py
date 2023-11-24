## Script to post-process the single-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter, alpha, from a data derived covariance matrix.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.

import numpy as np
import sys, os
from warnings import warn
from utils import symmetrized, cov_filter_smu, load_matrices_single, check_eigval_convergence, add_cov_terms, fit_shot_noise_rescaling, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from collect_raw_covariance_matrices import load_raw_covariances_smu


def compute_disconnected_term(EEaA1: np.ndarray[float], EEaA2: np.ndarray[float], RRaA1: np.ndarray[float], RRaA2: np.ndarray[float], RR: np.ndarray[float], weights: np.ndarray[float]):
    n_bins = len(RR)
    w_aA1 = RRaA1 / np.sum(RRaA1, axis=0)
    w_aA2 = RRaA2 / np.sum(RRaA2, axis=0)
    diff1 = EEaA1 - w_aA1 * EEaA1.sum(axis=0)
    diff2 = EEaA2 - w_aA2 * EEaA2.sum(axis=0)
    RRaRRb = np.matmul(np.asmatrix(RR).T, np.asmatrix(RR))
    fact = np.ones((n_bins, n_bins)) - np.matmul(np.asmatrix(weights).T, np.asmatrix(weights))
    cx = np.asarray(np.matmul(diff1.T, diff2) / np.matmul(fact, RRaRRb))
    return symmetrized(cx)

def optimize_alpha_jackknife(jackknife_file: str, weight_dir: str, input_file: dict, cov_filter: np.ndarray[int], tracer: int = 1, print_function = print) -> (float, int, np.ndarray[float], np.ndarray[float], np.ndarray[float], np.ndarray[float]):
    # Load jackknife xi estimates from data
    print_function(f"Loading correlation function jackknife estimates from {jackknife_file}")
    xi_jack = np.loadtxt(jackknife_file, skiprows=2)
    n_bins = xi_jack.shape[1] # total bins
    n_jack = xi_jack.shape[0] # total jackknives
    n = n_bins // m # radial bins

    weight_file = os.path.join(weight_dir, f'jackknife_weights_n{n}_m{m}_j{n_jack}_{tracer}{tracer}.dat')
    RR_file = os.path.join(weight_dir, f'binned_pair_counts_n{n}_m{m}_j{n_jack}_{tracer}{tracer}.dat')

    print_function("Loading weights file from %s"%weight_file)
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
    mean_xi = np.sum(xi_jack * weights, axis=0)
    tmp = weights * (xi_jack - mean_xi)
    data_cov = np.matmul(tmp.T, tmp)
    denom = np.matmul(weights.T, weights)
    data_cov /= (np.ones_like(denom) - denom)
    data_cov = data_cov[cov_filter]

    print_function(f"Loading weights file from {RR_file}")
    RR = np.loadtxt(RR_file)

    # Load in full jackknife theoretical matrices
    print_function("Loading best estimate of jackknife covariance matrix")
    c2j, c3j, c4j = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = True)
    c4j += compute_disconnected_term(input_file[f"EE1_{tracer}{tracer}_full"], input_file[f"EE2_{tracer}{tracer}_full"], input_file[f"RR1_{tracer}{tracer}_full"], input_file[f"RR2_{tracer}{tracer}_full"], RR, weights)[cov_filter]

    # Check matrix convergence
    check_eigval_convergence(c2j, c4j, "Jackknife")

    # Load in partial jackknife theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = True)
    c4s += np.array([compute_disconnected_term(EE1, EE2, RR1, RR2, RR, weights)[cov_filter]] for (EE1, EE2, RR1, RR2) in zip(input_file[f"EE1_{tracer}{tracer}"], input_file[f"EE2_{tracer}{tracer}"], input_file[f"RR1_{tracer}{tracer}"], input_file[f"RR2_{tracer}{tracer}"]))

    # Now optimize for shot-noise rescaling parameter alpha
    print_function("Optimizing for the shot-noise rescaling parameter")
    alpha_best = fit_shot_noise_rescaling(data_cov, c2j, c3j, c4j, c2s, c3s, c4s)
    print_function("Optimization complete - optimal rescaling parameter is %.6f"%alpha_best)

    # Compute jackknife covariance and precision matrices
    jack_cov = add_cov_terms(c2j, c3j, c4j, alpha_best)
    partial_jack_cov = add_cov_terms(c2s, c3s, c4s, alpha_best)
    _, jack_prec = compute_D_precision_matrix(partial_jack_cov, jack_cov)

    return alpha_best, n_jack, data_cov, jack_cov, partial_jack_cov, jack_prec

def post_process_jackknife(jackknife_file: str, weight_dir: str, file_root: str, n: int, m: int, outdir: str, skip_r_bins: int = 0, tracer: int = 1, print_function = print) -> dict[str]:
    cov_filter = cov_filter_smu(n, m, skip_r_bins)
    
    input_file = load_raw_covariances_smu(file_root, n, m, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    alpha_best, n_jack, data_cov, jack_cov, partial_jack_cov, jack_prec = optimize_alpha_jackknife(jackknife_file, weight_dir, input_file, cov_filter, tracer, print_function, full_return = True)

    # Load full covariance matrix terms
    c2f, c3f, c4f = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)
    # Compute full covariance matrix
    full_cov = add_cov_terms(c2f, c3f, c4f, alpha_best)

    # Check convergence
    check_eigval_convergence(c2f, c4f, "Full")

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial jackknife theoretical matrices
    c2fs, c3fs, c4fs = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = False)
    partial_cov = add_cov_terms(c2fs, c3fs, c4fs, alpha_best)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")    

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)

    output_dict = {"jackknife_theory_covariance": jack_cov, "full_theory_covariance": full_cov, "jackknife_data_covariance": data_cov, "shot_noise_rescaling": alpha_best, "jackknife_theory_precision": jack_prec, "full_theory_precision": full_prec, "N_eff": N_eff_D, "full_theory_D_matrix": full_D_est, "individual_theory_covariances": partial_cov, "individual_theory_jackknife_covariances": partial_jack_cov}

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Jackknife_n%d_m%d_j%d.npz' % (n, m, n_jack))
    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (6, 7):
        print("Usage: python post_process_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {OUTPUT_DIR} [{SKIP_R_BINS}]")
        sys.exit(1)

    jackknife_file = str(sys.argv[1])
    weight_dir = str(sys.argv[2])
    file_root = str(sys.argv[3])
    n = int(sys.argv[4])
    m = int(sys.argv[5])
    outdir = str(sys.argv[6])
    from utils import get_arg_safe
    skip_r_bins = get_arg_safe(7, int, 0)

    post_process_jackknife(jackknife_file, weight_dir, file_root, n, m, outdir, skip_r_bins)