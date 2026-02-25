"""
Functions to produce a combined/concatenated tracer covariances from two-tracer, as an alternative to the single-tracer approach used for e.g. `Valcin et al 2025 <https://arxiv.org/abs/2508.05467>`_.
Use with great care if the tracer combination pipeline involves weight rescaling and/or FKP weight updates.
"""

from .comb.convert_cov_legendre_multi_to_cat import convert_cov_legendre_multi_to_cat
from .comb.combine_covs_legendre_multi_to_cat import combine_covs_legendre_multi_to_cat