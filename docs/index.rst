RascalC Documentation
===================================

Overview 
----------

RascalC is a code primarily written in C++ to estimate covariance matrices from galaxy correlation functions. Given an input set of random particle locations and a correlation function, RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter. This is done by dividing the particles into jackknife regions, which are used to calibrate the rescaling parameter. RascalC can also be used to compute cross-covariances between different correlation functions.

The main estimators are described in O'Connell et al. 2016 (`arXiv <https://arxiv.org/abs/1510.01740>`_), O'Connell & Eisenstein 2018 (`arXiv <https://arxiv.org/abs/1808.05978>`_) and Philcox et al. 2019 (in prep.), with the final paper discussing this C++ implementation.

The source code is publicly available on `Github <https://github.com/oliverphilcox/RascalC>`_ and builds upon the Python package Rascal (`Github <https://github.com/rcoconnell/Rascal>`_).


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
   usage/visualization

For any queries regarding the code please contact `Oliver Philcox  <mailto:ohep2@alumni.cam.ac.uk>`_. 


.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
