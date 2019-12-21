RascalC: Fast Estimation of Galaxy Covariance Matrices
=======================================================

Overview
----------

RascalC is a code to quickly estimate covariance matrices from two- or three-point galaxy correlation functions, written in C++ and Python. Given an input set of random particle locations and a two-point correlation function (or input set of galaxy positions), RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter. For the 2PCF, the rescaling parameter can be calibrated by dividing the particles into jackknife regions and comparing sample to theoretical jackknife covariance. RascalC can also be used to compute Legendre-binned covariances and cross-covariances between different two-point correlation functions.

The main estimators are described in `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_, `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_ , `Philcox et al. 2019 <https://arxiv.org/abs/1904.11070>`_ and `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_ with the third and fourth papers discussing the new algorithms and C++ implementation. RascalC was also used in Philcox & Eisenstein (2019, accepted by MNRAS, `arXiv <https://arxiv.org/abs/1912.01010>`_) to compute the covariance of configuration-space power spectrum estimators.

The source code is publicly available on `Github <https://github.com/oliverphilcox/RascalC>`_ and builds upon the Python package `Rascal <https://github.com/rcoconnell/Rascal>`_. For general usage, comprehensive tutorials (:doc:`tutorial` and :doc:`tutorial_periodic`) are provided.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage/installation
   usage/getting-started
   usage/tutorial
   usage/tutorial_periodic
   usage/pre-processing
   usage/jackknife-weights
   usage/geometry-correction
   usage/correlation-functions
   usage/main-code
   usage/post-processing

For any queries regarding the code please contact `Oliver Philcox  <mailto:ohep2@alumni.cam.ac.uk>`_.


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
