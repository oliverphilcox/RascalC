Essential functions
===================

The main computation wrapper
----------------------------

This Python wrapper of ``RascalC`` is heavily interfaced with ``pycorr`` `library for 2-point correlation function estimation <https://github.com/cosmodesi/pycorr>`_.
Many of the arguments are intentionally similar to ``pycorr.TwoPointCorrelationFunction`` `high-level interface <https://py2pcf.readthedocs.io/en/latest/api/api.html#pycorr.correlation_function.TwoPointCorrelationFunction>`_.

Please bear with the long description; you can pay less attention to settings labeled optional in the beginning.

.. autofunction:: RascalC.run_cov


Post-processing
---------------

A suitable post-processing routine is invoked at the end of the main wrapper function (``RascalC.run_cov``), so in many circumstances you may not need to run it separately.
However, this automated but customizable post-processing routine is useful for timed-out runs, switching the mode, testing different cuts and/or output combinations in cases of insufficient convergence, etc.

.. autofunction:: RascalC.post_process_auto


Loading and exporting the final covariance matrices
---------------------------------------------------

.. automodule:: RascalC.cov_utils
    :members:
    :member-order: bysource