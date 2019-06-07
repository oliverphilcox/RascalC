RascalC: Fast Estimation of Galaxy Covariance Matrices
=======================================================

Overview 
----------

RascalC is a code to quickly estimate covariance matrices from galaxy 2-point correlation functions, written in C++ and Python. Given an input set of random particle locations and a correlation function (or input set of galaxy positions), RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter. This is done by dividing the particles into jackknife regions, which are used to calibrate the rescaling parameter. RascalC can also be used to compute cross-covariances between different correlation functions.

The main estimators are described in `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_, `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_ and `Philcox et al. 2019 <https://arxiv.org/abs/1904.11070>`_, with the final paper discussing the new algorithm and C++ implementation.

The source code is publicly available on `Github <https://github.com/oliverphilcox/RascalC>`_ and builds upon the Python package `Rascal <https://github.com/rcoconnell/Rascal>`_. For general usage, a comprehensive :doc:`tutorial` is provided.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   usage/installation
   usage/getting-started
   usage/tutorial
   usage/pre-processing
   usage/jackknife-weights
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
