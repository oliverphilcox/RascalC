.. RascalC documentation master file, created by
   sphinx-quickstart on Fri Jan 11 19:13:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

RascalC Documentation
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   pre-processing
   main-code
   utilities



RascalC is a code primarily written in C++ to estimate covariance matrices from galaxy correlation functions. Given an input set of random particle locations and a correlation function, RascalC produces an estimate of the associated covariance for a given binning strategy, with non-Gaussianities approximated by a 'shot-noise-rescaling' parameter. It can also be used to compute cross-covariances between different correlation functions.

The main estimators are described in O'Connell et al. 2016 ([arXiv](https://arxiv.org/abs/1510.01740)), O'Connell & Eisensein 2018 ([arXiv](https://arxiv.org/abs/1808.05978)) and Philcox et al. 2019 (in prep.), with the final paper discussing this implementation in C++.

The source code is publicly available on [Github](https://github.com/oliverphilcox/RascalC) and builds upon the Python package Rascal ([Github](https://github.com/rcoconnell/Rascal)).


Overview of RascalC
---------------------

Overview here

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
