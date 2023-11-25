# Contains some utility functions widely used in other scripts
# Not intended for execution from command line
import sys, os
import numpy as np
from warnings import warn
from astropy.io import fits
from scipy.optimize import fmin

def get_arg_safe(index: int, type = str, default: object = None) -> object:
    # get argument by index from sys.argv and convert it to the requested type if there are enough elements there
    # otherwise return the default value
    return type(sys.argv[index]) if len(sys.argv) > index else default

def blank_function(*args, **kwargs) -> None:
    # function that accepts anything and does nothing
    # mostly intended for skipping optional printing
    pass

def my_a2s(a, fmt='%.18e'):
    # custom array to string function
    return ' '.join([fmt % e for e in a])

def my_str_to_bool(s: str) -> bool:
    # naive conversion to bool for all non-empty strings is True, and one can't give an empty string as a command line argument, so need to make it more explicit
    return s not in ("0", "false")

def transposed(A: np.ndarray):
    # swap last two (matrix) axes
    return A.swapaxes(-2, -1)

def symmetrized(A: np.ndarray):
    # symmetrize a 2+D matrix over the last two axes
    return 0.5 * (A + transposed(A))

def parse_FKP_arg(FKP_weights: str) -> bool | tuple[float, str]:
    if not my_str_to_bool(FKP_weights): return False
    # determine if it actually has P0,NZ_name format. Such strings should convert to True
    arg_FKP_split = FKP_weights.split(",")
    if len(arg_FKP_split) == 2:
        return (float(arg_FKP_split[0]), arg_FKP_split[1])
    if len(arg_FKP_split) == 1: return True
    raise ValueError("FKP parameter matched neither USE_FKP_WEIGHTS (true/false in any register or 0/1) nor P0,NZ_name (float and string without space).")

def read_particles_fits_file(input_file: str, FKP_weights: bool | (float, str) = False, mask: int = 0, use_weights: bool = True):
    # Read FITS file with particles. Can apply mask filtering and compute FKP weights in different ways. Works for DESI setups
    filt = True # default pre-filter is true
    with fits.open(input_file) as f:
        data = f[1].data
        all_ra = data["RA"]
        all_dec = data["DEC"]
        all_z = data["Z"]
        colnames = data.columns.names
        all_w = data["WEIGHT"] if "WEIGHT" in colnames and use_weights else np.ones_like(all_z)
        if FKP_weights:
            all_w *= 1/(1+FKP_weights[0]*data[FKP_weights[1]]) if FKP_weights != True else data["WEIGHT_FKP"]
        if "WEIGHT" not in colnames and not FKP_weights: warn("No weights found, assigned unit weight to each particle.")
        if mask: filt = (data["STATUS"] & mask == mask) # all 1-bits from mask have to be set in STATUS; skip if mask=0
    return np.array((all_ra, all_dec, all_z, all_w)).T[filt]

def read_xi_file(xi_file: str):
    # Interpret RascalC text format using numpy functions
    if not os.path.isfile(xi_file): raise FileNotFoundError('Could not find input file %s' % xi_file)
    r_vals = np.genfromtxt(xi_file, max_rows=1)
    mu_vals = np.genfromtxt(xi_file, max_rows=1, skip_header=1)
    xi_vals = np.genfromtxt(xi_file, skip_header=2)
    return r_vals, mu_vals, xi_vals

def write_xi_file(xi_file: str, r_vals: np.ndarray[float], mu_vals: np.ndarray[float], xi_vals: np.ndarray[float]):
    # Reproduce RascalC text format using numpy functions
    header = my_a2s(r_vals) + '\n' + my_a2s(mu_vals)
    np.savetxt(xi_file, xi_vals, header=header, comments='')

def write_binning_file(out_file: str, r_edges: np.ndarray[float], print_function = blank_function):
    # Save bin edges array into a Corrfunc (and RascalC) radial binning file format
    np.savetxt(out_file, np.array((r_edges[:-1], r_edges[1:])).T)
    print_function("Binning file '%s' written successfully." % out_file)

def cov_filter_smu(n: int, m: int, skip_r_bins: int = 0):
    indices_1d = np.arange(m * skip_r_bins, m * n)
    return np.ix_(indices_1d, indices_1d)

def cov_filter_legendre(n: int, max_l: int, skip_r_bins: int = 0, skip_l: int = 0):
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    l_indices = np.arange(n_l - skip_l)
    r_indices = np.arange(skip_r_bins, n)
    indices_l_r = (n_l * r_indices)[:, None] + l_indices[None, :]
    # indices_l_r = (n_l * r_indices)[None, :] + l_indices[:, None] # could switch the ordering right here easily but then saved binary RascalC results will become incompatible
    indices_1d = indices_l_r.ravel()
    return np.ix_(indices_1d, indices_1d)

def load_matrices_single(input_data: dict[str], cov_filter: np.ndarray[int], tracer: int = 1, full: bool = True, jack: bool = False) -> (np.ndarray[float], np.ndarray[float], np.ndarray[float]):
    joint = "j" * jack + "_"
    suffix = + str(tracer) * 2 + "_full" * full
    c2 = input_data["c2" + joint + suffix]
    c3 = input_data["c3" + joint + str(tracer) + "," + suffix]
    c4 = input_data["c4" + joint + str(tracer) * 2 + "," + suffix]
    if full: # 2D matrices, filter can be applied directly
        c2 = symmetrized(c2[cov_filter])
        c3 = symmetrized(c3[cov_filter])
        c4 = symmetrized(c4[cov_filter])
    else: # 3D matrices, need to loop over the first index first
        c2 = np.array([symmetrized(_[cov_filter]) for _ in c2])
        c3 = np.array([symmetrized(_[cov_filter]) for _ in c3])
        c4 = np.array([symmetrized(_[cov_filter]) for _ in c4])
    return c2, c3, c4

def check_eigval_convergence(c2: np.ndarray[float], c4: np.ndarray[float], kind: str = "") -> None:
    eig_c4 = np.linalg.eigvalsh(c4)
    eig_c2 = np.linalg.eigvalsh(c2)
    if min(eig_c4) < -min(eig_c2):
        if kind and not kind.endswith(" "): kind += " "
        warn(f"{kind}4-point covariance matrix has not converged properly via the eigenvalue test. Min eigenvalue of C4 = {min(eig_c4):.2e}, min eigenvalue of C2 = {min(eig_c2):.2e}")

def check_positive_definiteness(full_cov: np.ndarray[float]) -> None:
    if np.any(np.linalg.eigvalsh(full_cov) <= 0): raise ValueError("The full covariance is not positive definite - insufficient convergence")

def add_cov_terms_single(c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], alpha: float = 1) -> np.ndarray[float]:
    return c4 + c3 * alpha + c2 * alpha**2

def compute_D_precision_matrix(partial_cov: np.ndarray[float], full_cov: np.ndarray[float]) -> (np.ndarray[float], np.ndarray[float]):
    n_samples = len(partial_cov)
    n_bins = len(full_cov)
    sum_partial_cov = np.sum(partial_cov, axis = 0)
    tmp = 0.
    for i in range(n_samples):
        c_excl_i = (sum_partial_cov - partial_cov[i]) / (n_samples - 1)
        tmp += np.matmul(np.linalg.inv(c_excl_i), partial_cov[i])
    full_D_est = (n_samples-1.) / n_samples * (-np.eye(n_bins) + tmp / n_samples)
    full_prec = np.matmul(np.eye(n_bins) - full_D_est, np.linalg.inv(full_cov))
    return full_D_est, full_prec

def compute_N_eff_D(full_D_est: np.ndarray[float], print_function = blank_function) -> float:
    n_bins = len(full_D_est)
    slogdetD = np.linalg.slogdet(full_D_est)
    D_value = slogdetD[0] * np.exp(slogdetD[1] / n_bins)
    if slogdetD[0] < 0:
        print_function("N_eff is negative! Setting to zero")
        N_eff_D = 0.
    else:
        N_eff_D = (n_bins+1.)/D_value+1.
        print_function("Total N_eff Estimate: %.4e"%N_eff_D)
    return N_eff_D

def Psi(alpha: float, c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], c2s: np.ndarray[float], c3s: np.ndarray[float], c4s: np.ndarray[float]):
    """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
    c_tot = add_cov_terms_single(c2, c3, c4, alpha)
    partial_covs = add_cov_terms_single(c2s, c3s, c4s, alpha)
    _, Psi = compute_D_precision_matrix(partial_covs, c_tot)
    return Psi

def neg_log_L1(alpha: float, target_cov: np.ndarray[float], c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], c2s: np.ndarray[float], c3s: np.ndarray[float], c4s: np.ndarray[float]):
    """Return negative log L1 likelihood between data and theory covariance matrices"""
    Psi_alpha = Psi(alpha, c2, c3, c4, c2s, c3s, c4s)
    logdet = np.linalg.slogdet(Psi_alpha)
    if logdet[0] < 0:
        # Remove any dodgy inversions
        return np.inf        
    return np.trace(np.matmul(Psi_alpha, target_cov)) - logdet[1]

def fit_shot_noise_rescaling(target_cov: np.ndarray[float], c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], c2s: np.ndarray[float], c3s: np.ndarray[float], c4s: np.ndarray[float]):
    alpha_best = fmin(neg_log_L1, 1., args = (target_cov, c2, c3, c4, c2s, c3s, c4s))
    return alpha_best

def load_matrices_multi(input_data: dict[str], cov_filter: np.ndarray[int], full: bool = True, jack: bool = False, ntracers: int = 2) -> (np.ndarray[float], np.ndarray[float], np.ndarray[float]):
    suffix_jack = "j" * jack
    suffix_full = "_full" * full
    single_array_shape = np.array(input_data["c4_11,11" + suffix_full].shape) # take reference shape from c4_11,11, with _full suffix if loading full
    single_array_shape[:-2] = np.array(np.zeros(single_array_shape[:-2])[cov_filter].shape) # account for the shape change with cov_filter
    single_array_shape = list(single_array_shape)

    # arrays to store cN with first N indices being tracers
    c2 = np.zeros([ntracers] * 2 + single_array_shape)
    c3 = np.zeros([ntracers] * 3 + single_array_shape)
    c4 = np.zeros([ntracers] * 4 + single_array_shape)
    # simpler arrays for accumulating number of additions to normalize by
    n2 = np.zeros([ntracers] * 2)
    n3 = np.zeros([ntracers] * 3)
    n4 = np.zeros([ntracers] * 4)

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
            t1, t2 = [int(c)-1 for c in index]
            # accumulate c2 with normal index order
            c2[t1, t2] += matrices
            n2[t1, t2] += 1
            # accumulate c2 with interchanged indices and transposed matrix, will ensure symmetry
            c2[t2, t1] += transposed(matrices)
            n2[t2, t1] += 1
        if matrix_name == "c3" + suffix_jack:
            if len(index) != 4: raise ValueError(f"Unexpected 3-point index length for {index}")
            if index[1] != ",": raise ValueError(f"Unexpected 3-point index format for {index}")
            t2, t1, t3 = int(index[0])-1, int(index[2])-1, int(index[3])-1
            # accumulate c3 with normal index order
            c3[t2, t1, t3] += matrices
            n3[t2, t1, t3] += 1
            # accumulate c3 with swapped j1 and j3
            c3[t2, t3, t1] += transposed(matrices)
            n3[t2, t3, t1] += 1
        if matrix_name == "c4" + suffix_jack:
            if len(index) != 5: raise ValueError(f"Unexpected 4-point index length for {index}")
            if index[2] != ",": raise ValueError(f"Unexpected 4-point index format for {index}")
            t1, t2, t3, t4 = int(index[0])-1, int(index[1])-1, int(index[3])-1, int(index[4])-1
            # All symmetries possible for c4 without transpositions
            permutations4 = ((t1, t2, t3, t4), # original
                            (t2, t1, t3, t4), # first two tracers interchanged
                            (t1, t2, t4, t3), # last two tracers interchanged
                            (t2, t1, t4, t3), # first and last two tracers interchanged at the same time
                            )
            for (i1, i2, i3, i4) in permutations4:
                c4[i1, i2, i3, i4] += matrices
                n4[i1, i2, i3, i4] += 1
                # now swap tracer pairs and transpose
                c4[i3, i4, i1, i2] += transposed(matrices)
                n4[i3, i4, i1, i2] += 1
        # else do nothing
    
    # check that all parts have been accumulated indeed
    if np.count_nonzero(n2 == 0) > 0: raise ValueError("Some 2-point terms missing")
    if np.count_nonzero(n3 == 0) > 0: raise ValueError("Some 3-point terms missing")
    if np.count_nonzero(n4 == 0) > 0: raise ValueError("Some 4-point terms missing")

    # now can normalize safely
    c2 /= n2
    c3 /= n4
    c4 /= n4

    return c2, c3, c4

def add_cov_terms_multi(c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], alphas: list[float] | np.ndarray[float], ntracers: int = 2) -> (np.ndarray[float], np.ndarray[float]):
    def construct_fields(t1: int, t2: int, t3: int, t4: int, alpha1: float, alpha2: float):
        # Reconstruct the full field for given input fields and rescaling parameters

        # Create kronecker deltas
        d_xw = (t1 == t4)
        d_xz = (t1 == t3)
        d_yw = (t2 == t4)
        d_yz = (t2 == t3)

        full = c4[t1, t2, t3, t4] + 0.25 * alpha1 * (d_xw * c3[t1, t2, t3] + d_xz * c3[t1, t2, t4]) + 0.25 * alpha2 * (d_yw * c3[t2, t1, t3] + d_yz * c3[t2, t1, t4]) + 0.5 * alpha1 * alpha2 * (d_xw * d_yz + d_xz * d_yw) * c2[t1, t2]
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

    single_array_shape = list(c2.shape)[2:]
    if c2.shape != [ntracers] * 2 + single_array_shape: raise ValueError("Unexpected shape of 2-point array")
    if c3.shape != [ntracers] * 3 + single_array_shape: raise ValueError("Unexpected shape of 3-point array")
    if c4.shape != [ntracers] * 4 + single_array_shape: raise ValueError("Unexpected shape of 4-point array")
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