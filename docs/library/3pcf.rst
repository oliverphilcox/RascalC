3PCF covariances (under development)
====================================

The main computation wrapper
----------------------------

This Python wrapper for ``RascalC`` 3PCF is interfaced with ``ENCORE`` <https://github.com/oliverphilcox/encore> for 3PCF ``pycorr`` <https://github.com/cosmodesi/pycorr> for 2PCF.
Many of the arguments are intentionally similar to ``pycorr.TwoPointCorrelationFunction`` `high-level interface <https://py2pcf.readthedocs.io/en/latest/api/api.html#pycorr.correlation_function.TwoPointCorrelationFunction>`_.

Please bear with the long description; you can pay less attention to settings labeled optional in the beginning.

.. autofunction:: RascalC.run_cov_3pcf


Post-processing
---------------

The post-processing routine is invoked automatically at the end of the main wrapper function (:func:`RascalC.run_cov_3pcf`), so in many circumstances you may not need to run it separately.
However, knowing the post-processing routine is useful for timed-out runs, switching the mode, testing different cuts and/or output combinations in cases of insufficient convergence, etc.

After post-processing, you probably want to run the extra convergence check (see :mod:`RascalC.convergence_check_extra`).

.. autofunction:: RascalC.post_process_3pcf.post_process_3pcf


Loading and exporting the final covariance matrices
---------------------------------------------------

.. automodule:: RascalC.cov_utils_3pcf
    :members:
    :member-order: bysource