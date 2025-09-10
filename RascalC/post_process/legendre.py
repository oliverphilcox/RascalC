## Script to post-process the single-field Legendre binned integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import os
from .utils import cov_filter_legendre, load_matrices_single, check_eigval_convergence, add_cov_terms_single, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_legendre
from typing import Literal, Callable, Iterable


def post_process_legendre(file_root: str, n: int, max_l: int, outdir: str, alpha: float = 1, skip_r_bins: int | tuple[int, int] = 0, skip_l: int = 0, tracer: Literal[1, 2] = 1, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print) -> dict[str]:
    cov_filter = cov_filter_legendre(n, max_l, skip_r_bins, skip_l)
    
    input_file = load_raw_covariances_legendre(file_root, n, max_l, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load in full theoretical matrices
    print_function("Loading best estimate of covariance matrix")
    c2, c3, c4 = load_matrices_single(input_file, cov_filter, tracer, full = True, jack = False)

    # Check matrix convergence
    check_eigval_convergence(c2, c4, alpha, print_function = print_function)

    # Compute full covariance matrices and precision
    full_cov = add_cov_terms_single(c2, c3, c4, alpha)

    # Check positive definiteness
    check_positive_definiteness(full_cov)

    # Compute full precision matrix
    print_function("Computing the full precision matrix estimate:")
    # Load in partial theoretical matrices
    c2s, c3s, c4s = load_matrices_single(input_file, cov_filter, tracer, full = False, jack = False)
    partial_cov = add_cov_terms_single(c2s, c3s, c4s, alpha)
    full_D_est, full_prec = compute_D_precision_matrix(partial_cov, full_cov)
    print_function("Full precision matrix estimate computed")

    # Now compute effective N:
    N_eff_D = compute_N_eff_D(full_D_est, print_function)

    output_dict = {"full_theory_covariance": full_cov, "shot_noise_rescaling": alpha, "full_theory_precision": full_prec, "N_eff": N_eff_D, "full_theory_D_matrix": full_D_est, "individual_theory_covariances": partial_cov}

    output_name = os.path.join(outdir, 'Rescaled_Covariance_Matrices_Legendre_n%d_l%d.npz'%(n,max_l))
    np.savez_compressed(output_name, **output_dict)
    output_dict["path"] = output_name
    output_dict["filename"] = os.path.basename(output_name)
    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict
