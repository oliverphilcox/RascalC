from pycorr import TwoPointCorrelationFunction
import lsstypes
import numpy as np
import numpy.typing as npt
from ..pycorr_utils.utils import reshape_pycorr
from ..lsstypes_utils.utils import reshape_lsstypes
from ..allcounts_utils import get_s_edges_from_allcounts, get_mu_edges_from_allcounts
from ..cov_utils import get_cov_header, load_cov_legendre_multi
from ..pycorr_utils.counts import get_counts_from_pycorr
from ..lsstypes_utils.counts import get_counts_from_lsstypes
from ..mu_bin_legendre_factors import compute_mu_bin_legendre_factors
from .utils import guess_allcounts_format
from typing import Callable, Literal


def combine_covs_legendre_multi_to_cat(rascalc_results1: str, rascalc_results2: str, allcounts_files1: list[str], allcounts_files2: list[str], output_cov_file: str, max_l: int, r_step: float = 1, skip_r_bins: int | tuple[int, int] = 0, bias1: float = 1, bias2: float = 1, output_cov_file1: str | None = None, output_cov_file2: str | None = None, allcounts_format: Literal[None, "pycorr", "lsstypes"] = None, print_function: Callable[[str], None] = print) -> npt.NDArray[np.float64]:
    """
    Given two-tracer Legendre mode RascalC results for two regions, produce a single-tracer covariance matrix for the combined/concatenated tracer (obtained by concatenating the catalogs of the two tracers, with weight in each optionally multiplied by the corresponding tracer's bias) for the region/footprint that is a combination of two regions/footprints.
    The correlations between the clustering statistics in the different regions are neglected, but the correlations between the two tracers in each region are included.
    More specifically, first the multi-tracer covariances for each region are converted to single-tracer covariances for the combined/concatenated tracer in each region, and then these single-tracer covariances are combined to produce the final covariance for the combined/concatenated tracer for the combined region.
    For additional details, see `Valcin et al 2025 <https://arxiv.org/abs/2508.05467>`_ and Appendix B.2 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.

    Parameters
    ----------
    rascalc_results1, rascalc_results2 : string
        Filenames for the RascalC (post-processing) results for the two regions in NumPy format.
    
    allcounts_files1, allcounts_files2 : list of strings
        Filenames for the ``pycorr`` (https://github.com/cosmodesi/pycorr) ``.npy`` or ``lsstypes`` (https://github.com/adematti/lsstypes) ``.h5``/``.hdf5``/``.txt`` files with the correlation functions and pair counts for the two regions.
        Each list must contain three filenames: first for the auto-correlation of the first tracer, second for the cross-correlation of the two tracers, and the third for the auto-correlation of the second tracer.
        The order of regions must be the same as in RascalC results.
    
    output_cov_file : string
        Filename for the output text file, in which the covariance matrix will be saved.

    max_l : integer
        The highest (even) multipole index, must match the RascalC results.

    r_step : float
        The width of the radial (separation) bins, must match the RascalC results.
    
    skip_r_bins : integer or tuple of two integers
        (Optional) removal of some radial bins from the loaded ``pycorr`` counts after adjusting the radial (separation) bin width to match the covariance settings.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.
        E.g. if the ``pycorr`` counts are in 1 Mpc/h bins from 0 to 200 Mpc/h and the RascalC covariances are computed only between 20 and 200 Mpc/h in 4 Mpc/h wide bins, ``skip_r_bins`` should be ``5`` or ``(5, 0)``.
    
    bias1, bias2 : float
        (Optional) the bias values to upweight the first and the second tracer respectively.
        Default is 1 for both tracers (i.e., no upweighting).
    
    output_cov_file1, output_cov_file2 : string or None
        (Optional) if provided, the text covariance matrices for the corresponding region for the combined/concatenated tracer will be saved in this file.

    allcounts_format : None, "pycorr" or "lsstypes"
        (Optional) the format of the allcounts files, either "pycorr" for files with pycorr TwoPointCorrelationFunction objects or "lsstypes" for files with lsstypes Count2Correlation objects. Default is None for auto-determination based on file extensions.
    
    print_function : Callable[[str], None]
        (Optional) custom function to use for printing. Needs to take string arguments and not return anything. Default is ``print``.

    Returns
    -------
    combined_cov : npt.NDArray[np.float64]
        The resulting covariance matrix for the combined/concatenated tracer in the combined region.
    """
    if len(allcounts_files1) != 3:
        raise ValueError("Need 3 allcounts_files1 for auto1, cross12, auto2")
    if len(allcounts_files2) != 3:
        raise ValueError("Need 3 allcounts_files2 for auto1, cross12, auto2")
    
    allcounts_format = guess_allcounts_format(allcounts_format, allcounts_files1 + allcounts_files2)
    
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov_legendre_multi(rascalc_results1, max_l, print_function)
    n_bins = len(cov1)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov_legendre_multi(rascalc_results2, max_l, print_function)

    # Read allcounts files to figure out weights
    allcounts1 = []
    for allcounts_file1 in allcounts_files1:
        if allcounts_format == "pycorr":
            allcounts1.append(reshape_pycorr(TwoPointCorrelationFunction.load(allcounts_file1), n_mu=None, r_step=r_step, skip_r_bins=skip_r_bins))
        else:
            allcounts1.append(reshape_lsstypes(lsstypes.read(allcounts_file1), n_mu=None, r_step=r_step, skip_r_bins=skip_r_bins))

    n_r_bins = len(get_s_edges_from_allcounts(allcounts1[0])) - 1
    mu_edges = get_mu_edges_from_allcounts(allcounts1[0])

    allcounts2 = []
    for allcounts_file2 in allcounts_files2:
        if allcounts_format == "pycorr":
            allcounts2.append(reshape_pycorr(TwoPointCorrelationFunction.load(allcounts_file2), n_mu=None, r_step=r_step, skip_r_bins=skip_r_bins))
        else:
            allcounts2.append(reshape_lsstypes(lsstypes.read(allcounts_file2), n_mu=None, r_step=r_step, skip_r_bins=skip_r_bins))

    # Add weighting by bias for each tracer
    bias_weights = (bias1**2, 2*bias1*bias2, bias2**2) # auto1, cross12, auto2 are multiplied by product of biases of tracers involved in each. Moreover, cross12 enters twice because wrapped cross21 is the same.
    if allcounts_format == "pycorr":
        # Now multiply ALL the counts
        allcounts1 = [allcounts * factor for allcounts, factor in zip(allcounts1, bias_weights)]
        allcounts2 = [allcounts * factor for allcounts, factor in zip(allcounts2, bias_weights)]
        # Extract not normalized RR counts
        weights1 = np.array([get_counts_from_pycorr(allcounts, counts_factor=1) for allcounts in allcounts1])
        weights2 = np.array([get_counts_from_pycorr(allcounts, counts_factor=1) for allcounts in allcounts2])
        # Other weights will be counts after summation and normalization
        weight1 = get_counts_from_pycorr(sum(allcounts1).normalize(), counts_factor=1)
        weight2 = get_counts_from_pycorr(sum(allcounts2).normalize(), counts_factor=1)
    else:
        # Extract not normalized RR counts and multiply by bias weights
        weights1 = np.array([get_counts_from_lsstypes(allcounts, counts_factor=1) * factor for allcounts, factor in zip(allcounts1, bias_weights)])
        weights2 = np.array([get_counts_from_lsstypes(allcounts, counts_factor=1) * factor for allcounts, factor in zip(allcounts2, bias_weights)])
        # Other weights will be counts after summation and normalization
        weight1 = np.sum(weights1, axis=0) / sum(allcounts.get('RR').norm * factor for (allcounts, factor) in zip(allcounts1, bias_weights)) * sum(allcounts.get(allcounts.count_names[0]).norm * factor for (allcounts, factor) in zip(allcounts1, bias_weights)) # first, take naive sum of not normalized counts multiplied by bias weights, then divide by the sum of the normalizations of RR counts multiplied by bias weights to get the correct normalization factor for the naively summed counts, and then multiply by the sum of the normalizations of the first counts (presumably DD) multiplied by bias weights. the latter mirrors the lsstypes code at https://github.com/adematti/lsstypes/blob/3bf32b393f81fa7068fbccd027fa793193e056c3/lsstypes/types.py#L1343
        weight2 = np.sum(weights2, axis=0) / sum(allcounts.get('RR').norm * factor for (allcounts, factor) in zip(allcounts2, bias_weights)) * sum(allcounts.get(allcounts.count_names[0]).norm * factor for (allcounts, factor) in zip(allcounts2, bias_weights))

    # Normalize weights
    sum_weight = weight1 + weight2
    weight1 /= sum_weight
    weight2 /= sum_weight
    # Normalize the full weights across correlation function labels
    weights1 /= np.sum(weights1, axis=0)[None, :, :]
    weights2 /= np.sum(weights2, axis=0)[None, :, :]

    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse=True)

    # First, convert multi-tracer cov to single-tracer in each region

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weights1/2[tracer, r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while tracers and radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,tkl,lj,km->tikjm', leg_mu_factors, weights1, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins // 3)
    pd2 = np.einsum('il,tkl,lj,km->tikjm', leg_mu_factors, weights2, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins // 3)
    # We have correct [t_in, l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape.
    # The resulting covs are single-tracer (for the combined catalogs) so there is no t_out.

    # Produce single-tracer covs for each region
    cov1 = pd1.T.dot(cov1).dot(pd1)
    cov2 = pd2.T.dot(cov2).dot(pd2)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header=header1) # includes shot-noise rescaling value in the header
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header=header2) # includes shot-noise rescaling value in the header
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    n_bins //= 3 # all covariances are single tracer now

    # Now, combine single-tracer covs

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weight1/2[r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd1 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight1, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins)
    pd2 = np.einsum('il,kl,lj,km->ikjm', leg_mu_factors, weight2, mu_leg_factors, np.eye(n_r_bins)).reshape(n_bins, n_bins)
    # We have correct [l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape

    # Produce and save combined cov
    cov = pd1.T.dot(cov1).dot(pd1) + pd2.T.dot(cov2).dot(pd2)
    np.savetxt(output_cov_file, cov, header=header) # includes source parts and their shot-noise rescaling values in the header
    return cov
