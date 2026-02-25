from pycorr import TwoPointCorrelationFunction
import numpy as np
from ..pycorr_utils.utils import reshape_pycorr
from ..cov_utils import get_cov_header, load_cov_legendre_multi
from ..pycorr_utils.counts import get_counts_from_pycorr
from ..mu_bin_legendre_factors import compute_mu_bin_legendre_factors
from typing import Callable


def load_cov_text(filename: str) -> tuple[np.typing.NDArray[np.float64], str]:
    cov = np.loadtxt(filename)
    # read header line if present
    header = '' # blank header by default
    with open(filename) as f:
        l = f.readline() # read the first line
        if l[0] == '#': # if it starts as a comment
            header = l[1:].strip() # take the rest of it as header, removing the leading/trailing spaces and the newline
    return cov, header

def convert_cov_legendre_multi_to_cat(rascalc_results: str, pycorr_files: list[str], output_cov_file: str, max_l: int, r_step: float = 1, skip_r_bins: int | tuple[int, int] = 0, bias1: float = 1, bias2: float = 1, print_function: Callable[[str], None] = print) -> np.typing.NDArray[np.float64]:
    """
    Given a two-tracer Legendre mode RascalC result (or a text covariance), produce a single-tracer covariance matrix for the combined/concatenated tracer (obtained by concatenating the catalogs of the two tracers, with weight in each optionally multiplied by the corresponding tracer's bias).
    The correlations between the two tracers in each region are included.
    For additional details, see `Valcin et al 2025 <https://arxiv.org/abs/2508.05467>`_ and Appendix B.2 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.

    Parameters
    ----------
    rascalc_results : string
        Filename for the RascalC two-tracer post-processing results in NumPy format, or the text file with the covariance matrix converted from such a NumPy file.
    
    pycorr_files : list of strings
        Filenames for the ``pycorr`` (https://github.com/cosmodesi/pycorr) ``.npy`` files with the correlation functions and pair counts.
        The list must contain three filenames: first for the auto-correlation of the first tracer, second for the cross-correlation of the two tracers, and the third for the auto-correlation of the second tracer.
    
    output_cov_file : string
        Filename for the output text file, in which the covariance matrix will be saved.

    max_l : integer
        The highest (even) multipole index, must match the RascalC results.

    r_step : float
        The width of the radial (separation) bins, must match the RascalC results.
    
    skip_r_bins : integer or tuple of two integers
        (Optional) removal of some radial bins from the loaded ``pycorr`` counts before adjusting the radial (separation) bin width to match the covariance settings.
        First (or the only) number sets the number of radial/separation bins to skip from the beginning.
        Second number (if provided) sets the number of radial/separation bins to skip from the end.
        By default, no bins are skipped.
        E.g. if the ``pycorr`` counts are in 1 Mpc/h bins from 0 to 200 Mpc/h and the RascalC covariances are computed only between 20 and 200 Mpc/h, ``skip_r_bins`` should be ``20``.
    
    bias1, bias2 : float
        (Optional) the bias values to upweight the first and the second tracer respectively.
        Default is 1 for both tracers (i.e., no upweighting).
    
    print_function : Callable[[str], None]
        (Optional) custom function to use for printing. Needs to take string arguments and not return anything. Default is ``print``.

    Returns
    -------
    combined_cov : np.typing.NDArray[np.float64]
        The resulting covariance matrix for the combined region.
    """
    # Read RascalC results
    if any(rascalc_results.endswith(ext) for ext in (".npy", ".npz")):
        # read numpy file
        header = get_cov_header(rascalc_results)
        cov_in = load_cov_legendre_multi(rascalc_results, max_l, print_function)
    else:
        # read text file
        header, cov_in = load_cov_text(rascalc_results)
    
    n_bins = len(cov_in)

    # Read pycorr files to figure out weights
    weights = []
    for pycorr_file in pycorr_files:
        xi_estimator = reshape_pycorr(TwoPointCorrelationFunction.load(pycorr_file), n_mu = None, r_step = r_step, skip_r_bins = skip_r_bins)
        weights.append(get_counts_from_pycorr(xi_estimator, counts_factor = 1))
    weights = np.array(weights)

    n = xi_estimator.shape[0]
    mu_edges = xi_estimator.edges[1]

    # Add weighting by bias for each tracer
    bias_weights = np.array((bias1**2, 2*bias1*bias2, bias2**2)) # auto1, cross12, auto2 are multiplied by product of biases of tracers involved in each. Moreover, cross12 enters twice because wrapped cross21 is the same.
    weights *= bias_weights[:, None, None]

    # Normalize weights across the correlation type axis
    weights /= np.sum(weights, axis=0)[None, :, :]

    mu_leg_factors, leg_mu_factors = compute_mu_bin_legendre_factors(mu_edges, max_l, do_inverse = True)

    # Derivatives of angularly binned 2PCF wrt Legendre are leg_mu_factors[ell//2, mu_bin]
    # Angularly binned 2PCF are added with weights (normalized) weights[tracer, r_bin, mu_bin]
    # Derivatives of Legendre wrt binned 2PCF are mu_leg_factors[mu_bin, ell//2]
    # So we need to sum such product over mu bins, while tracers and radial bins stay independent, and the partial derivative of combined 2PCF wrt the 2PCFs 1/2 will be
    pd = np.einsum('il,tkl,lj,km->tikjm', leg_mu_factors, weights, mu_leg_factors, np.eye(n)).reshape(n_bins, n_bins // 3)
    # We have correct [t_in, l_in, r_in, l_out, r_out] ordering and want to make these matrices in the end thus the reshape.
    # The output cov is single-tracer (for the combined catalog) so there is no t_out.

    # Produce and save combined cov
    cov_out = pd.T.dot(cov_in).dot(pd)
    np.savetxt(output_cov_file, cov_out, header = header)
    return cov_out
