import lsstypes
import numpy as np
from warnings import warn


def fix_bad_bins_lsstypes(xi_estimator: lsstypes.Count2Correlation) -> lsstypes.Count2Correlation:
    # fixes bins with negative wcounts by overwriting their content by reflection
    # only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32
    kw = {}
    for count_name in xi_estimator.count_names:
        counts : lsstypes.Count2 = xi_estimator.get(count_name)
        bad_bins_mask = counts.values('normalized_counts') < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {count_name}.normalized_counts ({counts.values('normalized_counts')[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.values('normalized_counts')[s_bin, -1-mu_bin]:.2e})")
            counts._data['normalized_counts'][s_bin, mu_bin] = counts._data['normalized_counts'][s_bin, -1-mu_bin]
            counts._data['norm'][s_bin, mu_bin] = counts._data['norm'][s_bin, -1-mu_bin]
        bad_bins_mask = counts.values('counts') < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {count_name}.counts ({counts.values('counts')[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.values('counts')[s_bin, -1-mu_bin]:.2e})")
            counts._data['normalized_counts'][s_bin, mu_bin] = counts._data['normalized_counts'][s_bin, -1-mu_bin]
            counts._data['norm'][s_bin, mu_bin] = counts._data['norm'][s_bin, -1-mu_bin]
        kw[count_name] = counts
    return lsstypes.Count2Correlation(**kw, estimator=xi_estimator.estimator, attrs=xi_estimator.attrs)
