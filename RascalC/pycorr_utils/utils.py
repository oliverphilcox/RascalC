import pycorr
import numpy as np
from warnings import warn
from ..xi.utils import write_xi_file # for convenience of use in other scripts


def fix_bad_bins_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator) -> pycorr.twopoint_estimator.BaseTwoPointEstimator:
    # fixes bins with negative wcounts by overwriting their content by reflection
    # only known cause for now is self-counts (DD, RR) in bin 0, n_mu_orig/2-1 – subtraction is sometimes not precise enough, especially with float32
    cls = xi_estimator.__class__
    kw = {}
    for name in xi_estimator.count_names:
        counts = getattr(xi_estimator, name)
        bad_bins_mask = counts.wcounts < 0
        for s_bin, mu_bin in zip(*np.nonzero(bad_bins_mask)):
            warn(f"Negative {name}.wcounts ({counts.wcounts[s_bin, mu_bin]:.2e}) found in bin {s_bin}, {mu_bin}; replacing them with reflected bin ({counts.wcounts[s_bin, -1-mu_bin]:.2e})")
            counts.wcounts[s_bin, mu_bin] = counts.wcounts[s_bin, -1-mu_bin]
        kw[name] = counts
    return cls(**kw)


def reshape_pycorr(xi_estimator: pycorr.twopoint_estimator.BaseTwoPointEstimator, n_mu: int | None = None, r_step: float | None = None, r_max: float = np.inf, skip_r_bins: int = 0) -> pycorr.twopoint_estimator.BaseTwoPointEstimator:
    n_mu_orig = xi_estimator.shape[1]
    if n_mu_orig % 2 != 0: raise ValueError("Wrapping not possible")
    if n_mu:
        if n_mu_orig % (2 * n_mu) != 0: raise ValueError("Angular rebinning not possible")
        mu_factor = n_mu_orig // 2 // n_mu
    else: mu_factor = 1 # leave the original number of mu bins

    if not r_step: r_factor = 1
    else:
        # determine the radius step in pycorr
        r_steps_orig = np.diff(xi_estimator.edges[0])
        r_step_orig = np.mean(r_steps_orig)
        if not np.allclose(r_steps_orig, r_step_orig, rtol=5e-3, atol=5e-3): raise ValueError("Radial rebinning only supported for linear bins")
        r_factor_exact = r_step / r_step_orig
        r_factor = int(np.rint(r_factor_exact))
        if not np.allclose(r_factor, r_factor_exact, rtol=5e-3): raise ValueError(f"Radial rebinning seems impossible: exact ratio of steps is {r_factor_exact}, closest integer is {r_factor} and that is too far")

    # Apply r_max cut
    r_values = xi_estimator.sepavg(axis = 0)
    r_bins_indices = np.where(r_values <= r_max)[0]
    xi_estimator = xi_estimator[:r_bins_indices.max() + 1] # can not apply mask to pycorr estimators; the bins are adjacent anyway

    return fix_bad_bins_pycorr(xi_estimator[skip_r_bins * r_factor:])[::r_factor, ::mu_factor].wrap() # first skip bins, then fix bad bins, then rebin and wrap to positive mu
