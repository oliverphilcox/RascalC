"This generates input xi for RascalC from lsstypes Count2Correlation objects/files"

import lsstypes
import numpy as np
import numpy.typing as npt


def get_input_xi_from_lsstypes(xi_estimator: lsstypes.Count2Correlation) -> npt.NDArray[np.float64]:
    # assume already wrapped; for input xi need to divide by SS instead of RR in post-recon case
    corr = xi_estimator.value()
    if 'SS' in xi_estimator.count_names: # for input xi need to divide by SS instead of RR in post-recon case, but SS may not be available in pre-recon case
        corr *= xi_estimator.get('RR').values('normalized_wcounts') / xi_estimator.get('SS').values('normalized_wcounts')
    return corr