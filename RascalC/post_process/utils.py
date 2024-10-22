import numpy as np
from warnings import warn
from ..utils import blank_function, symmetrized, transposed
from scipy.optimize import fmin
from scipy.linalg import sqrtm


def cov_filter_smu(n: int, m: int, skip_r_bins: int = 0):
    """Produce a 2D indexing array for s,mu covariance matrices."""
    indices_1d = np.arange(m * skip_r_bins, m * n)
    return np.ix_(indices_1d, indices_1d)


def cov_filter_legendre(n: int, max_l: int, skip_r_bins: int = 0, skip_l: int = 0):
    """Produce a 2D indexing array for Legendre covariance matrices."""
    if max_l % 2 != 0: raise ValueError("Only even multipoles supported")
    n_l = max_l // 2 + 1
    l_indices = np.arange(n_l - skip_l)
    r_indices = np.arange(skip_r_bins, n)
    indices_l_r = (n_l * r_indices)[:, None] + l_indices[None, :]
    # indices_l_r = (n_l * r_indices)[None, :] + l_indices[:, None] # could switch the ordering right here easily but then saved binary RascalC results will become incompatible
    indices_1d = indices_l_r.ravel()
    return np.ix_(indices_1d, indices_1d)


def load_matrices_single(input_data: dict[str], cov_filter: np.ndarray[int], tracer: int = 1, full: bool = True, jack: bool = False) -> tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]]:
    """Load the single-tracer covariance matrix terms. Allows to choose the tracer index (1-based) to load."""
    joint = "j" * jack + "_"
    suffix = str(tracer) * 2 + "_full" * full
    c2 = input_data["c2" + joint + suffix]
    c3 = input_data["c3" + joint + str(tracer) + "," + suffix]
    c4 = input_data["c4" + joint + str(tracer) * 2 + "," + suffix]
    matrices = (c2, c3, c4)

    def finalize_matrix(a: np.ndarray):
        return symmetrized(a[cov_filter])
    
    if full: # 2D matrices, filter can be applied directly
        c234 = [finalize_matrix(a) for a in matrices]
    else: # 3D matrices, need to loop over the first index first
        c234 = [np.array(list(map(finalize_matrix, a))) for a in matrices]
    return tuple(c234)


def check_eigval_convergence(c2: np.ndarray[float], c4: np.ndarray[float], alpha: float = 1, kind: str = "") -> bool:
    """Perform the eigenvalue convergence test on the covariance matrix terms.
    The default assumption is `alpha >= 1`.
    The condition is stronger for `alpha < 1`."""
    inv_sqrt_c2 = np.linalg.inv(sqrtm(c2))
    eig = np.linalg.eigvalsh(inv_sqrt_c2.dot(c4).dot(inv_sqrt_c2))
    if min(eig) <= -alpha**2:
        if kind and not kind.endswith(" "): kind += " "
        warn(f"{kind}4-point covariance matrix has not converged properly via the eigenvalue test for shot-noise rescaling >= {alpha}. Min eigenvalue of C2^{{-1/2}} C4 C2^{{-1/2}} = {min(eig):.2f}, should be > {-alpha**2:.2f}")
        return False
    return True


def check_positive_definiteness(full_cov: np.ndarray[float]) -> None:
    """Ensure that the final covariance matrix is positive definite; raise an error if it is not."""
    if np.any(np.linalg.eigvalsh(full_cov) <= 0): raise ValueError("The full covariance is not positive definite - insufficient convergence")


def add_cov_terms_single(c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], alpha: float = 1) -> np.ndarray[float]:
    """Add the single-tracer covariance matrix terms with a given shot-noise rescaling value."""
    return c4 + c3 * alpha + c2 * alpha**2


def compute_D_precision_matrix(partial_cov: np.ndarray[float], full_cov: np.ndarray[float]) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """Compute the quadratic order bias correction D and the precision bias with this correction."""
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
    """Compute the effective number of mocks (giving an equivalent covariance matrix precision) from the quadratic order bias correction factor D."""
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
    """Return negative log L1 likelihood between theory and target (data jackknife or mock sample) covariance matrices.
    log L1 is the Kullback-Leibler divergence with constant terms (including log(det(target_cov))) removed.
    As a result, the `target_cov` can be a singular matrix.
    This function does not allow negative shot-noise rescaling `alpha` by returning infinity."""
    if alpha < 0: return np.inf # negative shot-noise rescaling causes problems and does not make sense
    Psi_alpha = Psi(alpha, c2, c3, c4, c2s, c3s, c4s)
    logdet = np.linalg.slogdet(Psi_alpha)
    if logdet[0] < 0:
        # Remove any dodgy inversions
        return np.inf        
    return np.trace(np.matmul(Psi_alpha, target_cov)) - logdet[1]


def fit_shot_noise_rescaling(target_cov: np.ndarray[float], c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], c2s: np.ndarray[float], c3s: np.ndarray[float], c4s: np.ndarray[float]):
    """Fit the covariance matrix model to `target_cov` to find the optimal shot-noise rescaling.
    `target_cov` can be a singular matrix."""
    alpha_best = fmin(neg_log_L1, 1., args = (target_cov, c2, c3, c4, c2s, c3s, c4s))
    return alpha_best


def load_matrices_multi(input_data: dict[str], cov_filter: np.ndarray[int], full: bool = True, jack: bool = False, ntracers: int = 2) -> tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]]:
    """Load the multi-tracer covariance matrix terms."""
    suffix_jack = "j" * jack
    suffix_full = "_full" * full
    single_array_shape = np.array(input_data["c4_11,11" + suffix_full].shape) # take reference shape from c4_11,11, with _full suffix if loading full
    single_array_shape[-2:] = np.zeros(single_array_shape[-2:])[cov_filter].shape # account for the shape change with cov_filter
    single_array_shape = list(single_array_shape)

    # arrays to store cN with first N indices being tracers
    c2 = np.zeros([ntracers] * 2 + single_array_shape)
    c3 = np.zeros([ntracers] * 3 + single_array_shape)
    c4 = np.zeros([ntracers] * 4 + single_array_shape)
    # simpler arrays for accumulating number of additions to normalize by
    n2 = np.zeros([ntracers] * 2)
    n3 = np.zeros([ntracers] * 3)
    n4 = np.zeros([ntracers] * 4)

    def finalize_matrix(a: np.ndarray):
        return a[cov_filter]

    # accumulate the values
    for matrix_name, matrices in input_data.items():
        matrix_name_split = matrix_name.split("_")
        if len(matrix_name_split) != 2 + full: continue # should skip full if not loading full, and skip subsamples if loading full
        if full and matrix_name_split[-1] != "full": continue # double-check for safety
        if full: # 2D matrix, can apply cov filter directly
            matrices = finalize_matrix(matrices)
        else: # 3D matrix, should loops over the first index first
            matrices = np.array(list(map(finalize_matrix, matrices)))
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
    c2 /= n2.reshape([ntracers] * 2 + [1] * len(single_array_shape))
    c3 /= n3.reshape([ntracers] * 3 + [1] * len(single_array_shape))
    c4 /= n4.reshape([ntracers] * 4 + [1] * len(single_array_shape))

    return c2, c3, c4


def gen_corr_tracers(ntracers: int = 2) -> list[tuple[int, int]]:
    """Generate the pairs of tracer indices (zero-based) for all the correlation functions.
    Assumes that the X-Y cross-correlation is identical to Y-X.
    For ntracers = 1, gives `[(0, 0)]`.
    For ntracers = 2, gives `[(0, 0), (0, 1), (1, 1)]`, consistent with the old convention.
    When a tracer is added, the list begins in the same way."""
    corr_tracers = []
    for t2 in range(ntracers):
        for t1 in range(t2+1):
            corr_tracers.append((t1, t2))

    n_corr = ntracers * (ntracers + 1) // 2
    assert len(corr_tracers) == n_corr, "Mismatch in correlation functions counting or indices generation"
    return corr_tracers


def add_cov_terms_multi(c2: np.ndarray[float], c3: np.ndarray[float], c4: np.ndarray[float], alphas: list[float] | np.ndarray[float], ntracers: int = 2) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """Add the multi-tracer covariance matrix terms with given shot-noise rescaling values."""

    def construct_fields(t1: int, t2: int, t3: int, t4: int, alpha1: float, alpha2: float):
        # Reconstruct the full field for given input fields and rescaling parameters

        # Create kronecker deltas
        d_xw = (t1 == t4)
        d_xz = (t1 == t3)
        d_yw = (t2 == t4)
        d_yz = (t2 == t3)

        full = c4[t1, t2, t3, t4] + 0.25 * alpha1 * (d_xw * c3[t1, t2, t3] + d_xz * c3[t1, t2, t4]) + 0.25 * alpha2 * (d_yw * c3[t2, t1, t3] + d_yz * c3[t2, t1, t4]) + 0.5 * alpha1 * alpha2 * (d_xw * d_yz + d_xz * d_yw) * c2[t1, t2]
        return full

    single_array_shape = list(c2.shape)[2:]
    if c2.shape != tuple([ntracers] * 2 + single_array_shape): raise ValueError("Unexpected shape of 2-point array")
    if c3.shape != tuple([ntracers] * 3 + single_array_shape): raise ValueError("Unexpected shape of 3-point array")
    if c4.shape != tuple([ntracers] * 4 + single_array_shape): raise ValueError("Unexpected shape of 4-point array")
    n_bins = single_array_shape[-1]
    if single_array_shape[-2] != n_bins: raise ValueError("Covariance matrices are not square")
    samples_shape = single_array_shape[:-2]
    if len(samples_shape) > 1: raise ValueError("Multiple sample axes not implemented")

    corr_tracers = gen_corr_tracers(ntracers)
    n_corr = len(corr_tracers)

    c_tot = np.zeros([n_corr] * 2 + samples_shape + [n_bins] * 2) # array with each individual covariance accessible
    c_comb = np.zeros(samples_shape + [n_corr * n_bins] * 2) # full array suitable for inversion

    for i_corr1, (t1, t2) in enumerate(corr_tracers):
        alpha1, alpha2 = alphas[t1], alphas[t2]
        for i_corr2, (t3, t4) in enumerate(corr_tracers):
            tmp = construct_fields(t1, t2, t3, t4, alpha1, alpha2)
            c_tot[i_corr1, i_corr2] = tmp
            if samples_shape: # need extra sample axis
                c_comb[:, i_corr1*n_bins:(i_corr1+1)*n_bins, i_corr2*n_bins:(i_corr2+1)*n_bins] = tmp
            else:
                c_comb[i_corr1*n_bins:(i_corr1+1)*n_bins, i_corr2*n_bins:(i_corr2+1)*n_bins] = tmp

    return c_tot, symmetrized(c_comb) # add all remaining symmetries