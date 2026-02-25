"""
These functions follow a particular procedure for combination of DESI NGC and SGC 2-point correlation function measurements into GCcomb (previously N and S into NS).
For more information, see Appendix B of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.
"""

from .comb.combine_covs import combine_covs
from .comb.combine_covs_legendre import combine_covs_legendre
from .comb.combine_covs_multi import combine_covs_multi
from .comb.combine_covs_legendre_multi import combine_covs_legendre_multi