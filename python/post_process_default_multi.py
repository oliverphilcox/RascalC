## Script to post-process the multi-field integrals computed by the C++ code.
## We output the theoretical covariance matrices, (quadratic-bias corrected) precision matrices and the effective number of samples, N_eff.

import numpy as np
import sys, os
from utils import transposed, symmetrized, cov_filter_smu, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from collect_raw_covariance_matrices import load_raw_covariances_smu


def load_matrices_multi(input_data: dict[str], cov_filter: np.ndarray[int], full: bool = False, jack: bool = False, ntracers: int = 2):
    suffix_jack = "j" * jack
    suffix_full = "_full" * full
    single_array_shape = list(input_data["c4_11,11" + suffix_full][cov_filter].shape) # take reference shape from c4_11,11, with _full suffix if loading full

    # arrays to store cN with first N indices being tracers
    c2s = np.zeros([ntracers] * 2 + single_array_shape)
    c3s = np.zeros([ntracers] * 3 + single_array_shape)
    c4s = np.zeros([ntracers] * 4 + single_array_shape)
    # simpler arrays for accumulating number of additions to normalize by
    n2s = np.zeros([ntracers] * 2)
    n3s = np.zeros([ntracers] * 3)
    n4s = np.zeros([ntracers] * 4)

    # accumulate the values
    for matrix_name, matrices in input_data.values():
        matrix_name_split = matrix_name.split("_")
        if len(matrix_name_split) != 2 + full: continue # should skip full if not loading full, and skip subsamples if loading full
        if full and matrix_name_split[-1] != "full": continue # double-check for safety
        if full: # 2D matrix, can apply cov filter directly
            matrices = matrices[cov_filter]
        else: # 3D matrix, should loops over the first index first
            matrices = np.array(_[cov_filter] for _ in matrices)
        matrix_name = matrix_name_split[0]
        index = matrix_name_split[1]
        if matrix_name == "c2" + suffix_jack:
            if len(index) != 2: raise ValueError(f"Wrong 2-point index length for {index}")
            j1, j2 = [int(c)-1 for c in index]
            # accumulate c2 with normal index order
            c2s[j1, j2] += matrices
            n2s[j1, j2] += 1
            # accumulate c2 with interchanged indices and transposed matrix, will ensure symmetry
            c2s[j1, j2] += transposed(matrices)
            n2s[j1, j2] += 1
        if matrix_name == "c3" + suffix_jack:
            if len(index) != 4: raise ValueError(f"Unexpected 3-point index length for {index}")
            if index[1] != ",": raise ValueError(f"Unexpected 3-point index format for {index}")
            j2, j1, j3 = int(index[0])-1, int(index[2])-1, int(index[3])-1
            # accumulate c3 with normal index order
            c3s[j2, j1, j3] += matrices
            n3s[j2, j1, j3] += 1
            # accumulate c3 with swapped j1 and j3
            c3s[j2, j1, j3] += transposed(matrices)
            n3s[j2, j1, j3] += 1
        if matrix_name == "c4" + suffix_jack:
            if len(index) != 5: raise ValueError(f"Unexpected 4-point index length for {index}")
            if index[2] != ",": raise ValueError(f"Unexpected 4-point index format for {index}")
            j1, j2, j3, j4 = int(index[0])-1, int(index[1])-1, int(index[3])-1, int(index[4])-1
            # All symmetries possible for c4 without transpositions
            permutations4 = ((j1, j2, j3, j4), # original
                            (j2, j1, j3, j4), # first two indices interchanged
                            (j1, j2, j4, j3), # last two indices interchanged
                            (j2, j1, j4, j3), # first and last two indices interchanged at the same time
                            )
            for (i1, i2, i3, i4) in permutations4:
                c4s[i1, i2, i3, i4] += matrices
                n4s[i1, i2, i3, i4] += 1
                # now swap indices and transpose
                c4s[i3, i4, i1, i2] += transposed(matrices)
                n4s[i3, i4, i1, i2] += 1
        # else do nothing
    
    # check that all parts have been accumulated indeed
    if np.count_nonzero(n2s == 0) > 0: raise ValueError("Some 2-point terms missing")
    if np.count_nonzero(n3s == 0) > 0: raise ValueError("Some 3-point terms missing")
    if np.count_nonzero(n4s == 0) > 0: raise ValueError("Some 4-point terms missing")

    # now can normalize safely
    c2s /= n2s
    c3s /= n4s
    c4s /= n4s

    return c2s, c3s, c4s

def add_cov_terms_multi(c2s: np.ndarray[float], c3s: np.ndarray[float], c4s: np.ndarray[float], alphas: list[float] | np.ndarray[float], ntracers: int = 2):
    def construct_fields(t1, t2, t3, t4, alpha1, alpha2):
        # Reconstruct the full field for given input fields and rescaling parameters

        # Create kronecker deltas
        d_xw = (t1 == t4)
        d_xz = (t1 == t3)
        d_yw = (t2 == t4)
        d_yz = (t2 == t3)

        full = c4s[t1, t2, t3, t4] + 0.25 * alpha1 * (d_xw * c3s[t1, t2, t3] + d_xz * c3s[t1, t2, t4]) + 0.25 * alpha2 * (d_yw * c3s[t2, t1, t3] + d_yz * c3s[t2, t1, t4]) + 0.5 * alpha1 * alpha2 * (d_xw * d_yz + d_xz * d_yw) * c2s[t1, t2]
        return full

    corr_tracers = []
    for t2 in range(ntracers):
        for t1 in range(t2+1):
            corr_tracers.append([t1, t2])
    # for ntracers = 1, gives [[0, 0]]
    # for ntracers = 2, gives [[0, 0], [0, 1], [1, 1]], consistent with the old convention
    # when a tracer is added, the list begins in the same way

    n_corr = ntracers * (ntracers + 1) // 2
    assert len(corr_tracers) == n_corr, "Mismatch in correlation functions counting or indices generation"

    single_array_shape = list(c2s.shape)[2:]
    if c2s.shape != [ntracers] * 2 + single_array_shape: raise ValueError("Unexpected shape of 2-point array")
    if c3s.shape != [ntracers] * 3 + single_array_shape: raise ValueError("Unexpected shape of 3-point array")
    if c4s.shape != [ntracers] * 4 + single_array_shape: raise ValueError("Unexpected shape of 4-point array")
    n_bins = single_array_shape[-1]
    if single_array_shape[-2] != n_bins: raise ValueError("Covariance matrices are not square")
    samples_shape = single_array_shape[:-2]
    if len(samples_shape) > 1: raise ValueError("Multiple sample axes not implemented")

    c_tot = np.zeros([n_corr] * 2 + samples_shape + [n_bins] * 2) # array with each individual covariance accessible
    c_comb = np.zeros(samples_shape + [n_corr * n_bins] * 2) # full array suitable for inversion

    for i_corr1 in range(n_corr):
        t1, t2 = corr_tracers[i_corr1]
        alpha1, alpha2 = alphas[t1], alphas[t2]
        for i_corr2 in range(n_corr):
            t3, t4 = corr_tracers[i_corr2]
            tmp = construct_fields(t1, t2, t3, t4, alpha1, alpha2)
            c_tot[i_corr1, i_corr2] = tmp
            if samples_shape: # need extra sample axis
                c_comb[:, i_corr1*n_bins:(i_corr1+1)*n_bins, i_corr2*n_bins:(i_corr2+1)*n_bins] = tmp
            else:
                c_comb[i_corr1*n_bins:(i_corr1+1)*n_bins, i_corr2*n_bins:(i_corr2+1)*n_bins] = tmp

    return c_tot, symmetrized(c_comb) # add all remaining symmetries

def post_process_default_multi(file_root: str, n: int, m: int, outdir: str, alpha_1: float = 1, alpha_2: float = 1, skip_r_bins: int = 0, print_function = print):
    cov_filter = cov_filter_smu(n, m, skip_r_bins)

    input_file = load_raw_covariances_smu(file_root, n, m, print_function)

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

    output_name =os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Default_n%d_m%d.npz'%(n,m))
    np.savez(output_name, **output_dict)

    print_function("Saved output covariance matrices as %s"%output_name)

    return output_dict

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (5, 7, 8):
        print("Usage: python post_process_default_multi.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2} [{SKIP_R_BINS}]]")
        sys.exit(1)

    file_root = str(sys.argv[1])
    n = int(sys.argv[2])
    m = int(sys.argv[3])
    outdir = str(sys.argv[4])
    from utils import get_arg_safe
    alpha_1 = get_arg_safe(5, float, 1)
    alpha_2 = get_arg_safe(6, float, 1)
    skip_r_bins = get_arg_safe(7, int, 0)

    post_process_default_multi(file_root, n, m, outdir, alpha_1, alpha_2, skip_r_bins)
