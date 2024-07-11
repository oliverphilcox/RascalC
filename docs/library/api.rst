API
=================

The main wrapper function
-------------------------

This Python wrapper of ``RascalC`` is heavily interfaced with ``pycorr`` `library for 2-point correlation function estimation <https://github.com/cosmodesi/pycorr>`_.
Many of the arguments are intentionally similar to ``pycorr.TwoPointCorrelationFunction`` `high-level interface <https://py2pcf.readthedocs.io/en/latest/api/api.html#pycorr.correlation_function.TwoPointCorrelationFunction>`_.

Please bear with the long description; you can pay less attention to settings labeled optional in the beginning.

.. autofunction:: RascalC.run_cov