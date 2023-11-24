## Script to post-process the single-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys, os
from utils import cov_filter_smu, load_matrices_single, check_eigval_convergence, add_cov_terms, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from collect_raw_covariance_matrices import load_raw_covariances_smu


def post_process_default(file_root: str, n: int, m: int, outdir: str, alpha: float = 1, skip_r_bins: int = 0, tracer: int = 1, print_function = print) -> dict[str]:
    cov_filter = cov_filter_smu(n, m, skip_r_bins)

    input_file = load_raw_covariances_smu(file_root, n, m, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c2, c3, c4 = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)

    # Check matrix convergence
    check_eigval_convergence(c2, c4)

    # Compute full covariance matrices and precision
    full_cov = add_cov_terms(c2, c3, c4, alpha)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = False)

    partial_cov = add_cov_terms(c2s, c3s, c4s, alpha)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)

    output_dict = {"full_theory_covariance": full_cov, "shot_noise_rescaling": alpha, "full_theory_precision": full_prec, "N_eff": N_eff_D, "full_theory_D_matrix": full_D_est, "individual_theory_covariances": partial_cov}

    output_name = os.path.join(outdir, f'Rescaled_Covariance_Matrices_Default_n{n}_m{m}.npz')
    np.savez_compressed(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (5, 6, 7):
        print("Usage: python post_process_default.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]")
        sys.exit(1)

    file_root = str(sys.argv[1])
    n = int(sys.argv[2])
    m = int(sys.argv[3])
    outdir = str(sys.argv[4])
    from utils import get_arg_safe
    alpha = get_arg_safe(5, float, 1)
    skip_r_bins = get_arg_safe(6, int, 0)
    
    post_process_default(file_root, n, m, outdir, alpha, skip_r_bins)