## Script to post-process the multi-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_smu, load_matrices_multi, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_smu
from typing import Iterable, Callable


def post_process_default_multi(file_root: str, n: int, m: int, outdir: str, alpha_1: float = 1, alpha_2: float = 1, skip_r_bins: int | tuple[int, int] = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print) -> dict[str]:
    cov_filter = cov_filter_smu(n, m, skip_r_bins)

    input_file = load_raw_covariances_smu(file_root, n, m, n_samples, print_function)

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

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Default_n%d_m%d.npz'%(n,m))
    np.savez_compressed(output_name, **output_dict)
    output_dict["path"] = output_name
    output_dict["filename"] = os.path.basename(output_name)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict
