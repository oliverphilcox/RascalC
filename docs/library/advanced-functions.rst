Advanced functions
==================

This list is not complete, but aims to provide references for the more used/useful extended features.


Covariance comparison functions
-------------------------------

.. automodule:: RascalC.cov_comparison
    :members:
    :member-order: bysource


Extra convergence check functions
---------------------------------

Extra convergence check is performed in the main wrapper function (``RascalC.run_cov``) and also performed by default in the automatic post-processing function (``RascalC.post_process_auto``), but the following functions may be useful to run additionally.

.. automodule:: RascalC.convergence_check_extra
    :members:
    :member-order: bysource


Sample catenation function
--------------------------

.. autofunction:: RascalC.cat_raw_covariance_matrices


Covariance combination functions for two independent regions
------------------------------------------------------------

These functions follow the particular procedure for combination of DESI NGC and SGC 2-point correlation function measurements into GCcomb (previously N and S into NS).
For more information, see Appendix B of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.

.. automodule:: RascalC.combine_regions
    :members:
    :imported-members:
    :member-order: bysource

Sample covariance utility functions
-----------------------------------

Radial and angular bins (``s_mu`` binning mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: RascalC.pycorr_utils.sample_cov
    :members:
    :member-order: bysource

Legendre multipoles (``legendre`` binning modes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: RascalC.pycorr_utils.sample_cov_multipoles
    :members:
    :member-order: bysource