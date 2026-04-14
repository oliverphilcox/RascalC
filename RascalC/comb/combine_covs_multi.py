from pycorr import TwoPointCorrelationFunction
import lsstypes
import numpy as np
import numpy.typing as npt
from ..pycorr_utils.utils import reshape_pycorr
from ..lsstypes_utils.utils import reshape_lsstypes
from ..cov_utils import get_cov_header, load_cov
from ..pycorr_utils.counts import get_counts_from_pycorr
from ..lsstypes_utils.counts import get_counts_from_lsstypes
from .utils import guess_allcounts_format
from typing import Callable, Literal


def combine_covs_multi(rascalc_results1: str, rascalc_results2: str, allcounts_files1: list[str], allcounts_files2: list[str], output_cov_file: str, n_mu_bins: int | None = None, r_step: float = 1, skip_r_bins: int | tuple[int, int] = 0, output_cov_file1: str | None = None, output_cov_file2: str | None = None, allcounts_format: Literal[None, "pycorr", "lsstypes"] = None, print_function: Callable[[str], None] = print) -> npt.NDArray[np.float64]:
    """
    Produce s,mu mode two-tracer covariance matrix for the region/footprint that is a combination of two regions/footprints neglecting the correlations between the clustering statistics in the different regions.
    For additional details, see Appendix B.1 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.

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

    n_mu_bins : integer
        The number of angular (mu) bins, must match the RascalC results.

    r_step : float
        The width of the radial (separation) bins, must match the RascalC results.
    
    skip_r_bins : integer or tuple of two integers
        (Optional) removal of some radial bins from the loaded ``pycorr`` counts after adjusting the radial (separation) bin width to match the covariance settings.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.
        E.g. if the ``pycorr`` counts are in 1 Mpc/h bins from 0 to 200 Mpc/h and the RascalC covariances are computed only between 20 and 200 Mpc/h in 4 Mpc/h wide bins, ``skip_r_bins`` should be ``5`` or ``(5, 0)``.
    
    output_cov_file1, output_cov_file2 : string or None
        (Optional) if provided, the text covariance matrices for the corresponding region will be saved in this file.

    allcounts_format : None, "pycorr" or "lsstypes"
        (Optional) the format of the allcounts files, either "pycorr" for files with pycorr TwoPointCorrelationFunction objects or "lsstypes" for files with lsstypes Count2Correlation objects. Default is None for auto-determination based on file extensions.

    print_function : Callable[[str], None]
        (Optional) custom function to use for printing. Needs to take string arguments and not return anything. Default is ``print``.

    Returns
    -------
    combined_cov : npt.NDArray[np.float64]
        The resulting covariance matrix for the combined region.
    """
    if len(allcounts_files1) != len(allcounts_files2): raise ValueError("Need the same number of allcounts files for both results")
    # Read RascalC results
    header1 = get_cov_header(rascalc_results1)
    cov1 = load_cov(rascalc_results1, print_function)
    header2 = get_cov_header(rascalc_results2)
    cov2 = load_cov(rascalc_results2, print_function)
    # Save to their files if any
    if output_cov_file1: np.savetxt(output_cov_file1, cov1, header=header1)
    if output_cov_file2: np.savetxt(output_cov_file2, cov2, header=header2)
    header = f"combined from {rascalc_results1} with {header1} and {rascalc_results2} with {header2}" # form the final header to include both

    allcounts_format = guess_allcounts_format(allcounts_format, allcounts_files1 + allcounts_files2)

    # Read allcounts files to figure out weights
    weight1 = np.zeros(0)
    for allcounts_file1 in allcounts_files1:
        if allcounts_format == "pycorr":
            weight1 = np.append(weight1, get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(allcounts_file1).normalize(), n_mu=n_mu_bins, r_step=r_step, skip_r_bins=skip_r_bins), counts_factor=1).ravel())
        else:
            xi_estimator = reshape_lsstypes(lsstypes.read(allcounts_file1), n_mu=n_mu_bins, r_step=r_step, skip_r_bins=skip_r_bins)
            weight1 = np.append(weight1, (get_counts_from_lsstypes(xi_estimator) * xi_estimator.get(xi_estimator.count_names[0]).norm).ravel()) # the first counts are presumably DD - mirroring the lsstypes code at https://github.com/adematti/lsstypes/blob/3bf32b393f81fa7068fbccd027fa793193e056c3/lsstypes/types.py#L1343
    weight2 = np.zeros(0)
    for allcounts_file2 in allcounts_files2:
        if allcounts_format == "pycorr":
            weight2 = np.append(weight2, get_counts_from_pycorr(reshape_pycorr(TwoPointCorrelationFunction.load(allcounts_file2).normalize(), n_mu=n_mu_bins, r_step=r_step, skip_r_bins=skip_r_bins), counts_factor=1).ravel())
        else:
            xi_estimator = reshape_lsstypes(lsstypes.read(allcounts_file2), n_mu=n_mu_bins, r_step=r_step, skip_r_bins=skip_r_bins)
            weight2 = np.append(weight2, (get_counts_from_lsstypes(xi_estimator) * xi_estimator.get(xi_estimator.count_names[0]).norm).ravel())

    # Produce and save combined cov
    # following xi = (xi1 * weight1 + xi2 * weight2) / (weight1 + weight2)
    cov = (cov1 * weight1[None, :] * weight1[:, None] + cov2 * weight2[None, :] * weight2[:, None]) / (weight1 + weight2)[None, :] / (weight1 + weight2)[:, None]
    np.savetxt(output_cov_file, cov, header=header) # includes source parts and their shot-noise rescaling values in the header
    return cov
