"""
Function to post-process the multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameters, alpha_i, from data derived covariance matrices.
We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. The effective number of samples, N_eff, is also computed.
"""

import numpy as np
import os
from warnings import warn
from .utils import gen_corr_tracers, cov_filter_smu, load_matrices_multi, check_eigval_convergence, fit_shot_noise_rescaling, add_cov_terms_multi, check_positive_definiteness, compute_D_precision_matrix, compute_N_eff_D
from ..raw_covariance_matrices import load_raw_covariances_smu
from typing import Iterable, Callable


def load_disconnected_term_multi(input_data: dict[str], cov_filter: np.typing.NDArray[np.int_], RR: np.typing.NDArray[np.float64] | list[np.typing.NDArray[np.float64]], weights: np.typing.NDArray[np.float64] | list[np.typing.NDArray[np.float64]], full: bool = True, ntracers: int = 2) -> tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]:
    suffix_full = "_full" * full

    disconnected_array_names = ["EE1", "RR1", "EE2", "RR2"]

    EE_RR_single_array_shape = list(input_data[disconnected_array_names[0] + "_11" + suffix_full][cov_filter].shape)

    disconnected_arrays = np.zeros([len(disconnected_array_names)] + [ntracers] * 2 + EE_RR_single_array_shape)

    # put EE1/2 and RR1/2 into arrays with tracer indices
    for matrix_name, matrices in input_data.items():
        matrix_name_split = matrix_name.split("_")
        if len(matrix_name_split) != 2 + full: continue # should skip full if not loading full, and skip subsamples if loading full
        if full and matrix_name_split[-1] != "full": continue # double-check for safety
        matrix_name = matrix_name_split[0]
        if matrix_name not in disconnected_array_names: continue # only process the disconnected arrays
        matrix_index = disconnected_array_names.index(matrix_name)
        index = matrix_name_split[1]
        if len(index) != 2: raise ValueError(f"Wrong 2-point index length for {index}")
        t1, t2 = [int(c)-1 for c in index]
        # assign the array with normal tracer order
        disconnected_arrays[matrix_index, t1, t2] = matrices
        # assign the array with swapped tracers, only makes sense if they are not the same
        if t1 != t2:
            disconnected_arrays[matrix_index, t2, t1] = matrices

    RR_full = np.zeros([ntracers] * 2 + list(RR[0].shape))
    weights_full = np.zeros([ntracers] * 2 + list(weights[0].shape))

    # put full RR and weights into arrays with tracer indices
    for i_corr, (t1, t2) in gen_corr_tracers(ntracers):
        # assign the arrays with normal tracer order
        RR_full[t1, t2] = RR[i_corr]
        weights_full[t1, t2] = weights[i_corr]
        # assign the arrays with swapped tracers, only makes sense if they are not the same
        if t1 != t2:
            RR_full[t2, t1] = RR[i_corr]
            weights_full[t2, t1] = weights[i_corr]

    c4_single_array_shape = np.array(input_data["c4_11,11" + suffix_full].shape) # take reference shape from c4_11,11, with _full suffix if loading full
    c4_single_array_shape[:-2] = np.array(np.zeros(c4_single_array_shape[:-2])[cov_filter].shape) # account for the shape change with cov_filter
    c4_single_array_shape = list(c4_single_array_shape)

    cx = np.zeros([ntracers] * 4 + c4_single_array_shape)

    # Finally, fill cx
    for t1 in range(ntracers):
        for t2 in range(ntracers):
            for t3 in range(ntracers):
                for t4 in range(ntracers):
                    RRaRRb = np.matmul(np.asmatrix(RR_full[t1, t2]).T, np.asmatrix(RR_full[t3, t4]))
                    fact = 1 - np.matmul(np.asmatrix(weights_full[t1, t2]).T, weights_full[t3, t4])
                    norm = RRaRRb * fact

                    def compute_disconnected_term(EEaA1: np.typing.NDArray[np.float64], RRaA1: np.typing.NDArray[np.float64], EEaA2: np.typing.NDArray[np.float64], RRaA2: np.typing.NDArray[np.float64]):
                        w_aA1 = RRaA1 / RRaA1.sum(axis = 0)
                        w_aA2 = RRaA2 / RRaA2.sum(axis = 0)
                        diff1 = EEaA1 - w_aA1 * EEaA1.sum(axis = 0)
                        diff2 = EEaA2 - w_aA2 * EEaA2.sum(axis = 0)
                        cx = np.matmul(diff1.T, diff2) / norm
                        return cx[cov_filter]
                    
                    def get_disconnected_terms(disconnected_arrays: np.typing.NDArray[np.float64]):
                        # applies symmetry
                        # uses the fact that first two arrays are EE1 and RR1 and the last two are EE2 and RR2
                        return 0.5 * (compute_disconnected_term(*disconnected_arrays[:2, t1, t2], *disconnected_arrays[-2:, t3, t4]) + compute_disconnected_term(*disconnected_arrays[-2:, t1, t2], *disconnected_arrays[:2, t3, t4]))
                        # note that tracer indices are NOT generally swappable because of the way the norm is fixed beforehand, while EE/RR 1/2 with the same indices are in principle interchangeable
                    
                    if full: # 2D arrays
                        cx[t1, t2, t3, t4] = get_disconnected_terms(disconnected_arrays)
                    else:
                        cx[t1, t2, t3, t4] = np.array(list(map(get_disconnected_terms, np.moveaxis(disconnected_arrays, 3, 0))))
                        # first axis of disconnected_arrays is the array type, the next two are tracer indices, then comes the subsample axis, which we put in front and loop over

    return cx


def post_process_jackknife_multi(jackknife_file_11: str, jackknife_file_12: str, jackknife_file_22: str, weight_dir: str, file_root: str, m: int, outdir: str, skip_r_bins: int | tuple[int, int] = 0, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function: Callable[[str], None] = print, dry_run: bool = False):
    ## First load jackknife xi estimates from data:
    print_function("Loading correlation function jackknife estimates")
    xi_jack_11 = np.loadtxt(jackknife_file_11, skiprows=2)
    xi_jack_12 = np.loadtxt(jackknife_file_12, skiprows=2)
    xi_jack_22 = np.loadtxt(jackknife_file_22, skiprows=2)
    if not (xi_jack_11.shape == xi_jack_22.shape == xi_jack_12.shape): raise ValueError('Must have the same configuration of jackknives for each field.')

    n_bins = xi_jack_11.shape[1] # total bins
    n_jack = xi_jack_11.shape[0] # total jackknives
    n = n_bins // m # radial bins

    output_name = os.path.join(outdir, 'Rescaled_Multi_Field_Covariance_Matrices_Jackknife_n%d_m%d_j%d.npz' % (n, m, n_jack))
    name_dict = dict(path=output_name, filename=os.path.basename(output_name))
    if dry_run: return name_dict

    # First exclude any dodgy jackknife regions
    good_jk = np.where(np.all(np.isfinite(xi_jack_11) & np.isfinite(xi_jack_12) & np.isfinite(xi_jack_22), axis=1))[0] # all xi in jackknife have to be normal numbers
    if len(good_jk) < n_jack:
        warn("Using only %d out of %d jackknives - some xi values were not finite" % (len(good_jk), n_jack))

    xi_all = np.zeros([len(good_jk),3*n_bins])
    weights_all = np.zeros_like(xi_all)

    # Load in all xi functions
    xi_all[:, :n_bins] = xi_jack_11[good_jk]
    xi_all[:, n_bins:2*n_bins] = xi_jack_12[good_jk]
    xi_all[:, 2*n_bins:] = xi_jack_22[good_jk]

    # Load in all weights:
    weight_file11 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack))
    weight_file12 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_12.dat'%(n,m,n_jack))
    weight_file22 = os.path.join(weight_dir, 'jackknife_weights_n%d_m%d_j%d_22.dat'%(n,m,n_jack))
    weights11 = np.loadtxt(weight_file11)[:, 1:]
    weights12 = np.loadtxt(weight_file12)[:, 1:]
    weights22 = np.loadtxt(weight_file22)[:, 1:]
    weights = np.array([these_weights[good_jk] for these_weights in (weights11, weights12, weights22)])
    weights /= np.sum(weights, axis = 1)[:, None, :] # renormalize after possibly discarding some jackknives
    weights_all = weights.swapaxes(0, 1).reshape(len(good_jk), 3*n_bins)

    # Compute full covariance matrix:
    tmp_cov = np.zeros([len(good_jk), 3*n_bins])
    mean_xi = np.sum(xi_all * weights_all, axis=0)
    tmp_cov = weights_all * (xi_all-mean_xi)

    print_function("Computing full data covariance matrix")
    # Now compute covariance matrix:
    num = np.matmul(tmp_cov.T, tmp_cov)
    denom = np.matmul(weights_all.T, weights_all)
    data_cov = num/(np.ones_like(denom)-denom)

    # Load in all RR counts:
    RR_file11 = os.path.join(weight_dir, 'binned_pair_counts_n%d_m%d_j%d_11.dat'%(n,m,n_jack))
    RR_file12 = os.path.join(weight_dir, 'binned_pair_counts_n%d_m%d_j%d_12.dat'%(n,m,n_jack))
    RR_file22 = os.path.join(weight_dir, 'binned_pair_counts_n%d_m%d_j%d_22.dat'%(n,m,n_jack))
    RR11 = np.loadtxt(RR_file11)[:, 1:]
    RR12 = np.loadtxt(RR_file12)[:, 1:]
    RR22 = np.loadtxt(RR_file22)[:, 1:]
    RRs = np.array((RR11, RR12, RR22))

    cov_filter = cov_filter_smu(n, m, skip_r_bins)

    input_file = load_raw_covariances_smu(file_root, n, m, n_samples, print_function)

    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load autocovariance from data-covariance
    data_cov_11 = data_cov[:n_bins,:n_bins]
    data_cov_22 = data_cov[2*n_bins:,2*n_bins:]
    auto_data_cov = [data_cov_11,data_cov_22]

    alpha_best = np.ones(2) # fill with ones by default, although this should not matter

    # Load full jack matrices
    c2j, c3j, c4j = load_matrices_multi(input_file, cov_filter, full = True, jack = True)
    c4j += load_disconnected_term_multi(input_file, cov_filter, RRs, weights, full = True)
    # Load subsample jack matrices
    c2s, c3s, c4s = load_matrices_multi(input_file, cov_filter, full = False, jack = True)
    c4s += load_disconnected_term_multi(input_file, cov_filter, RRs, weights, full = False)

    ## Optimize for alpha_1 and alpha_2 separately.
    for t, this_data_cov in enumerate(auto_data_cov):

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

    np.savez_compressed(output_name, **output_dict)
    print_function("Saved output covariance matrices as %s" % output_name)

    return output_dict | name_dict
