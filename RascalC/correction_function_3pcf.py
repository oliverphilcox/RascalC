import os
import numpy as np
import numpy.typing as npt
from typing import Callable
from scipy.special import legendre
from .correction_function import compute_V_n_w_bar, load_randoms


def compute_inv_phi_periodic_3pcf(n: int, n_multipoles: int) -> npt.NDArray[np.float64]:
    "Compute the inverse 3PCF survey correction function multipoles for the periodic box geometry."
    ## Output periodic survey correction function
    phi_inv_mult = np.zeros([n, n, n_multipoles])
    
    ## Set to correct periodic survey values
    phi_inv_mult[:, :, 0] = 1

    return phi_inv_mult


def check_triple_counts_positive(leg_triple: npt.NDArray[np.float64], lenient_samebins: bool = False, print_function: Callable[[str], None] = print) -> None:
    "Check for negative counts, which should be problematic"
    n_mu = 2001
    triple_counts = np.zeros(list(leg_triple.shape[:-1]) + [n_mu])
    mu_values = np.linspace(-1, 1, n_mu)
    for ell in range(leg_triple.shape[-1]):
        triple_counts += leg_triple[:, :, ell][:, :, None] * legendre(ell)(mu_values)[None, None, :]
    problem_indices = np.argwhere(triple_counts <= 0)
    if len(problem_indices) > 0:
        rbins, mu_counts = np.unique(problem_indices[:, :2], axis=0, return_counts=True)
        for ((rbin1, rbin2), mu_count) in zip(rbins, mu_counts):
            if rbin1 > rbin2: continue # this case can be skipped by symmetry
            title = "WARNING"
            if lenient_samebins and rbin1 == rbin2: title = "INFO" # the problem for same-bin pairs is less critical (and seems more likely), for different-bin pairs it is more critical
            print_function(f"{title}: counts are not positive for radial bin pair {rbin1}, {rbin2} for {mu_count} mu values of {n_mu} checked")


def check_inv_phi_values(phi_inv_mult: npt.NDArray[np.float64], print_function: Callable[[str], None] = print) -> None:
    "Check that the mean of the monopole is neither too small nor too large"
    print_function(f"INFO: mean±std of the monopole of the inverse survey correction function is {np.mean(phi_inv_mult[:, :, 0])}±{np.std(phi_inv_mult[:, :, 0], ddof=1)}, expected to be not very far from 1. NB: this diagnostic depends on an empirical volume estimate with ConvexHull, which may be off (overestimated) by a factor of a few for realistic survey geometries.")
    if np.mean(phi_inv_mult[:, :, 0]) < 1e-3:
        raise ValueError("Survey correction function seems too small - are the RRR counts normalized correctly?")
    if np.mean(phi_inv_mult[:, :, 0]) > 1e3:
        raise ValueError("Survey correction function seems too large - are the RRR counts normalized correctly?")


def compute_inv_phi_aperiodic_3pcf(n: int, m: int, n_multipoles: int, r_bins: npt.NDArray[np.float64], triple_counts: npt.NDArray[np.float64], print_function: Callable[[str], None] = print) -> npt.NDArray[np.float64]:
    "Compute the inverse 3PCF survey correction function multipoles for the realistic survey geometry."

    mu_all = np.linspace(-1,1,m+1)
    mu_cen = 0.5*(mu_all[1:]+mu_all[:-1])
    
    ## reshape RRR counts and add symmetries
    RRR_true = triple_counts.reshape(n, n, m)
    RRR_true = (RRR_true + RRR_true.transpose(1, 0, 2)) / 2 # although wouldn't they be symmetric in radial bins already, having come from triple_counts?
        
    ## Now construct Legendre moments
    leg_triple = np.zeros([n, n, n_multipoles])
    for ell in range(n_multipoles):
        # (NB: we've absorbed a factor of delta_mu into RRR_true here)
        leg_triple[:, :, ell] += (2.*ell+1.) * np.sum(legendre(ell)(mu_cen)[None, None, :] * RRR_true, axis=-1) # shouldn't this be divided by 2 due to the legendre polynomial normalization?
    
    # as a precaution, check for negative counts
    check_triple_counts_positive(leg_triple, n_multipoles, print_function=print_function)

    vol_r = 4 * np.pi / 3 * (r_bins[:, 1] **3 - r_bins[:, 0] ** 3)

    ## Construct multipoles of inverse Phi
    phi_inv_mult = leg_triple / (.5 * vol_r[:, None, None] * vol_r[None, :, None]) # wouldn't it be easier to remove 0.5 from here and use 3 instead of 6 in the normalization?
            
    ## Check all seems reasonable
    check_inv_phi_values(phi_inv_mult, print_function=print_function)

    return phi_inv_mult


def compute_3pcf_correction_function(randoms_pos: npt.NDArray[np.float64], randoms_weights: npt.NDArray[np.float64], binfile: str, outdir: str, periodic: bool, RRR_file: str | None = None, print_function: Callable[[str], None] = print) -> str:
    """
    Function to compute the multipole decomposition of the 3PCF inverse survey correction function.
    The 3PCF survey correction function is defined as the ratio between idealistic and true RRR pair counts for a single survey.

    NB: Input RRR counts should be normalized by the cube of the sum of random weights here.
    NB: Assume mu is in [-1,1] limit here
    """

    n_multipoles = 7 # matches the value hard-coded in the C++ code

    if periodic:
        print_function("Assuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere")
    elif RRR_file is None:
        raise TypeError("RRR counts file is required if aperiodic")

    V, n_bar, w_bar = compute_V_n_w_bar(randoms_pos, randoms_weights)

    # Load in binning files 
    r_bins = np.loadtxt(binfile)
    n=len(r_bins)

    ## Define normalization constant
    norm = 6. * V * n_bar**3 * w_bar**3 # I don't think there is an exactly right answer once number density or weights vary across the survey
    # this is equal to 6 * np.sum(randoms_weights)**3 / V**2, but the other form is easier to compare against the definitions in the paper (Sections 4.2 and 5.2.1 of https://arxiv.org/pdf/1910.04764)

    print_function("Normalizing output survey correction by %.2e"%norm)

    if periodic:
        phi_inv_mult = compute_inv_phi_periodic_3pcf(n, n_multipoles)

    else:
        triple_counts = np.loadtxt(RRR_file)*np.sum(randoms_weights)**3
    
        # Compute number of angular bins in data-set
        m = (len(triple_counts)//n)//n
        if len(triple_counts) % m != 0: raise ValueError("Incorrect RRR format")

        phi_inv_mult = compute_inv_phi_aperiodic_3pcf(n, m, n_multipoles, r_bins, triple_counts / norm, print_function=print_function)
        
    if periodic:
        outfile = os.path.join(outdir, 'BinCorrectionFactor3PCF_n%d_periodic.txt'%(n))
    else:
        outfile = os.path.join(outdir, 'BinCorrectionFactor3PCF_n%d_m%d.txt'%(n,m))
    
    np.savetxt(outfile, phi_inv_mult.reshape(n*n, n_multipoles) * norm)
    print_function("Saved (normalized) output to %s\n"%outfile)

    return outfile


def compute_3pcf_correction_function_from_files(random_filename: str, binfile: str, outdir: str, periodic: bool, RRR_file: str | None = None, print_function: Callable[[str], None] = print) -> str:
    print_function("Loading randoms")
    return compute_3pcf_correction_function(*load_randoms(random_filename), binfile, outdir, periodic, RRR_file, print_function = print_function)


def compute_3pcf_correction_function_from_encore(randoms_pos: npt.NDArray[np.float64], randoms_weights: npt.NDArray[np.float64], binfile: str, outdir: str, triple_counts: npt.NDArray[np.float64], print_function: Callable[[str], None] = print) -> str:
    """
    Function to compute the multipole decomposition of the 3PCF inverse survey correction function from ENCORE triple counts.
    The 3PCF survey correction function is defined as the ratio between idealistic and true RRR pair counts for a single survey.

    NB: Input RRR counts are not normalized here, and are already in multipole format.
    Caveat: ENCORE only computes the triple counts for pairs of different radial bins, whereas RascalC expects them for all pairs of radial bins, crucially including the pairs of identical bins. Here we try to fill the missing data for those identical-bin pairs using the neighboring bin pairs. These should only affect the covariance rows and columns corresponding to the identical-bin pairs, which should be removed from the covariance for use with ENCORE 3PCF measurements in the end. So, those missing bin pairs should not matter in the end, but it is nicer not to have complete nonsence in the intermediate products.
    """

    n_multipoles = 7 # matches the value hard-coded in the C++ code

    ells = np.arange(len(triple_counts)) # rows in triple_counts correspond to multipoles, columns - to radial bin pairs; the first column might just contain these ells

    if np.array_equal(triple_counts[:, 0], ells): triple_counts = triple_counts[:, 1:] # remove the first column if it is all ells

    # Load in binning files 
    r_bins = np.loadtxt(binfile)
    n = len(r_bins)

    if triple_counts.shape[1] != n*(n-1)//2: raise ValueError("The shape of RRR_counts is inconsistent with the radial bins provided")
    # check this after removing the ells column if present
    # the columns correspond to radial bins

    triple_counts = 6 * triple_counts # the factor of 6 comes from the fact that ENCORE counts each triplet only once, whereas RascalC counts each triplet 6 times (once for each permutation of the three points). both should encounter every triplet in every relevant permutations on top of this. don't use *= to avoid overriding the original argument outside the function
    # change normalization from ENCORE to simple Legendre polynomials used in RascalC (according to Equation 4.11 in https://arxiv.org/pdf/1910.04764)
    triple_counts *= ((-1)**ells * np.sqrt(2 * ells + 1) / (4 * np.pi))[:, None] # add the second dimension, corresponding to the radial bins, to avoid indexing errors. should be fine to use *= now because triple_counts is a new local variable
    # the ell-dependent factor between the ENCORE 3-point basis functions and Legendre polynomials given by Equation (17) in https://arxiv.org/pdf/2105.08722
    # need to check if it is not division; there might also be a factor of 2 or something similar

    # ensure the number of multipoles in triple_counts is right to avoid indexing errors
    if len(triple_counts) < n_multipoles: # this seems more likely
        print_function(f"INFO: ENCORE triple counts have {len(triple_counts)} multipoles, fewer than {n_multipoles} used for the survey correction function, extending by zeros")
        triple_counts = np.vstack([triple_counts, np.zeros([n_multipoles - len(triple_counts), triple_counts.shape[1]])])
    elif len(triple_counts) > n_multipoles: # this seems less likely
        print_function(f"INFO: ENCORE triple counts have {len(triple_counts)} multipoles, more than {n_multipoles} used for the survey correction function, discarding the higher multipoles")
        triple_counts = triple_counts[:n_multipoles]

    bin_indices = np.arange(n)
    bin_index1 = np.repeat(bin_indices, n-1-bin_indices)
    bin_index2 = np.concatenate([bin_indices[i+1:] for i in range(n)])
    # bin_index1 and bin_index2 cover all the bin pairs under the condition bin_index1 < bin_index2, the order follows the ENCORE format
    # they could be read from first two non-comment rows of the ENCORE file, but this seems unnecessary

    leg_triple = np.zeros([n, n, n_multipoles])
    leg_triple[bin_index1, bin_index2] = triple_counts.T # fill above the diagonal; transposition puts radial bin pair index first and multipole index last
    leg_triple[bin_index2, bin_index1] = triple_counts.T # fill below the diagonal symmetrically
    # leave zeros at the diagonal for now

    vol_r = 4 * np.pi / 3 * (r_bins[:, 1] ** 3 - r_bins[:, 0] ** 3) # volume of radial/separation bins as 1D array

    V, n_bar, w_bar = compute_V_n_w_bar(randoms_pos, randoms_weights)

    ## Define normalization constant
    norm = 6. * V * n_bar**3 * w_bar**3 # I don't think there is an exactly right answer once number density or weights vary across the survey
    # this is equal to 6 * np.sum(randoms_weights)**3 / V**2, but the other form is easier to compare against the definitions in the paper (Sections 4.2 and 5.2.1 of https://arxiv.org/pdf/1910.04764)

    ## Construct multipoles of inverse Phi
    phi_inv_mult = leg_triple / (.5 * norm * vol_r[:, None, None] * vol_r[None, :, None])

    # fill the middle diagonal elements, which have been zeros
    # seems better to do in phi_inv_mult because its values change less
    bin_indices_middle = bin_indices[1:-1]
    phi_inv_mult[bin_indices_middle, bin_indices_middle] = (phi_inv_mult[bin_indices_middle+1, bin_indices_middle] + phi_inv_mult[bin_indices_middle-1, bin_indices_middle]) / 2 # average the neighboring elements along the column. the neighboring elements along the row are the same due to symmetry

    # fill the edge/corner diagonal elements, which are more tricky and have been zeros
    phi_inv_mult[0, 0] = (2 * (2 * phi_inv_mult[1, 0] - phi_inv_mult[2, 0]) + (2 * phi_inv_mult[1, 1] - phi_inv_mult[2, 2])) / 3
    phi_inv_mult[-1, -1] = (2 * (2 * phi_inv_mult[-2, -1] - phi_inv_mult[-3, -1]) + (2 * phi_inv_mult[-2, -2] - phi_inv_mult[-3, -3])) / 3
    
    # check for negative triple counts (or rather the correction function as it is passed to RascalC), which should be problematic
    check_triple_counts_positive(phi_inv_mult * norm, lenient_samebins=True, print_function=print_function)
            
    ## Check all seems reasonable
    check_inv_phi_values(phi_inv_mult, print_function=print_function)
        
    outfile = os.path.join(outdir, 'BinCorrectionFactor3PCF_n%d.txt' % n)
    
    np.savetxt(outfile, phi_inv_mult.reshape(n*n, n_multipoles) * norm)
    print_function("Saved (normalized) output to %s\n"%outfile)

    return outfile