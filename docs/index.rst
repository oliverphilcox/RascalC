RascalC: Fast Estimation of Galaxy Covariance Matrices
=======================================================

Overview
----------

RascalC is a code to quickly estimate covariance matrices from two- or three-point galaxy correlation functions, written in C++ and Python. Given an input set of random particle locations and a two-point correlation function (or input set of galaxy positions), RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter. For the 2PCF, the rescaling parameter can be calibrated by dividing the particles into jackknife regions and comparing sample to theoretical jackknife covariance. RascalC can also be used to compute Legendre-binned covariances and cross-covariances between different two-point correlation functions.

The main estimators are described in `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_, `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_ , `Philcox et al. 2019 <https://arxiv.org/abs/1904.11070>`_ and `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_ with the third and fourth papers discussing the new algorithms and C++ implementation. `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ is dedicated to applications to `DESI <https://desi.lbl.gov>`_ and updated validation techniques. RascalC was also used in Philcox & Eisenstein (2019, accepted by MNRAS, `arXiv <https://arxiv.org/abs/1912.01010>`_) to compute the covariance of configuration-space power spectrum estimators.

The source code is publicly available on `Github <https://github.com/oliverphilcox/RascalC>`_ and builds upon the Python package `Rascal <https://github.com/rcoconnell/Rascal>`_. For general usage, comprehensive tutorials (:doc:`legacy/tutorial` and :doc:`legacy/tutorial_periodic`) are provided.

.. toctree::
   :maxdepth: 2
   :caption: Python library

   library/installation
   library/essential-functions
   library/tutorials-examples
   library/advanced-functions

.. toctree::
   :maxdepth: 2
   :caption: Legacy (command line) usage

   legacy/installation
   legacy/getting-started
   legacy/tutorial
   legacy/tutorial_periodic
   legacy/pre-processing
   legacy/jackknife-weights
   legacy/geometry-correction
   legacy/correlation-functions
   legacy/main-code
   legacy/post-processing

For any queries regarding the code please contact `Michael 'Misha' Rashkovetskyi  <mailto:mrashkovetskyi@cfa.harvard.edu>`_.


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
