RascalC: Fast Estimation of Galaxy Covariance Matrices
=======================================================

Overview
----------

RascalC is a code to quickly estimate covariance matrices for two- or three-point galaxy correlation functions, written in C++ and Python.
Given an input set of random particle locations and a two-point correlation function (or input set of galaxy positions), RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter.
For the 2PCF, the rescaling parameter can be calibrated by dividing the particles into jackknife regions and comparing sample to theoretical jackknife covariance.
RascalC can also be used to compute Legendre-binned covariances and cross-covariances between different two-point correlation functions.

The main estimators are described in `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_, `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_, `Philcox et al. 2020 <https://arxiv.org/abs/1904.11070>`_ and `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_ with the third and fourth papers discussing the new algorithms and C++ implementation.
`Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_ is dedicated to applications to `DESI <https://desi.lbl.gov>`_ and updated validation techniques.
`Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_ introduces a new mode for Legendre-binned covariances and presents an extended validation using DESI DR1 mocks.
These two papers are combined (with repetitions eliminated) in Chapter 2 of `Rashkovetskyi 2025 dissertation <https://rashkovetsky.im/files/Dissertation_Rashkovetskyi_2025.pdf>`_.
RascalC was also used in `Philcox & Eisenstein 2020 <https://arxiv.org/abs/1912.01010>`_ to compute the covariance of configuration-space power spectrum estimators.

The source code is publicly available on `Github <https://github.com/oliverphilcox/RascalC>`_ and builds upon the Python package `Rascal <https://github.com/rcoconnell/Rascal>`_.
RascalC is now also available as a Python package (see :doc:`library/installation`).
For general usage, we provide :doc:`library/essential-functions` and :doc:`library/tutorials-examples`; you can also refer to legacy-mode/command-line :doc:`legacy/tutorial` and :doc:`legacy/tutorial_periodic`.

.. toctree::
   :maxdepth: 2
   :caption: Python library (2PCF)

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
