## Script to post-process the multi-field Legendre binned integrals computed by the C++ code, given a shot-noise rescaling parameter alpha.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys,os
from utils import cov_filter_legendre, load_matrices_multi, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from collect_raw_covariance_matrices import load_raw_covariances_legendre


def post_process_legendre_multi(file_root: str, n: int, max_l: int, outdir: str, alpha_1: float = 1, alpha_2: float = 1, skip_r_bins: int = 0, skip_l: int = 0, print_function = print):
    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)

    input_file = load_raw_covariances_legendre(file_root, n, max_l, print_function)

    alphas = [alpha_1, alpha_2]

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load full matrices
    c2, c3, c4 = load_matrices_multi(input_file, cov_filter, full = True, jack = False)
    c_tot, c_comb = add_cov_terms_multi(c2, c3, c4, alphas)

    # Check positive definiteness
    check_positive_definiteness(c_comb)

    # Load subsampled matrices (all submatrices combined)
    c2s, c3s, c4s = load_matrices_multi(input_file, cov_filter, full = False, jack = False)
    _, c_comb_subsamples = add_cov_terms_multi(c2s, c3s, c4s, alphas)

    # Now compute precision matrix
    D_est, prec_comb = compute_D_precision_matrix(c_comb_subsamples, c_comb)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff = compute_N_eff_D(D_est, print_function)

    output_dict = {"full_theory_covariance": c_comb, "all_covariances": c_tot, "shot_noise_rescaling": alphas, "full_theory_precision": prec_comb, "N_eff": N_eff, "full_theory_D_matrix": D_est, "individual_theory_covariances": c_comb_subsamples}

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Legendre_n%d_l%d.npz'%(n,max_l))

    np.savez(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (5, 7, 8, 9):
        print("Usage: python post_process_legendre_multi.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2}]")
        sys.exit(1)

    file_root = str(sys.argv[1])
    n = int(sys.argv[2])
    max_l = int(sys.argv[3])
    outdir = str(sys.argv[4])
    from utils import get_arg_safe
    alpha_1 = get_arg_safe(5, float, 1)
    alpha_2 = get_arg_safe(6, float, 1)
    skip_r_bins = get_arg_safe(7, int, 0)
    skip_l = get_arg_safe(8, int, 0)

    post_process_legendre_multi(file_root, n, max_l, outdir, alpha_1, alpha_2, skip_r_bins, skip_l)
