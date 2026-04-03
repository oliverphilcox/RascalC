"This generates binned pair counts for RascalC from lsstypes Count2Correlation objects/files"

import lsstypes
import numpy as np
import numpy.typing as npt


def get_counts_from_lsstypes(xi_estimator: lsstypes.Count2Correlation, counts_factor: float | None = None, split_above: float = np.inf) -> npt.NDArray[np.float64]:
    if not counts_factor: # use normalized counts
        return xi_estimator.get('RR').values('normalized_counts')
    paircounts = xi_estimator.get('R1R2').values('counts') / counts_factor
    nonsplit_mask = (xi_estimator.coords('s') < split_above)
    if split_above > 0: paircounts[nonsplit_mask] /= counts_factor # divide once more below the splitting scale
    return paircounts


def convert_counts_from_lsstypes_to_file(xi_estimator: lsstypes.Count2Correlation, outfile_name: str, counts_factor: float | None = None, split_above: float = np.inf) -> None:
    # xi_estimator_reshaped = reshape_lsstypes(xi_estimator, n_mu, r_step, r_max)
    paircounts = get_counts_from_lsstypes(xi_estimator, counts_factor, split_above)
    ## Write to file using numpy funs
    np.savetxt(outfile_name, paircounts.reshape(-1, 1)) # the file always has 1 column


def convert_counts_from_lsstypes_files(infile_name: str, outfile_name: str, counts_factor: float | None = None, split_above: float = np.inf) -> None:
    xi_estimator_orig = lsstypes.read(infile_name)
    convert_counts_from_lsstypes_to_file(xi_estimator_orig, outfile_name, counts_factor, split_above)
