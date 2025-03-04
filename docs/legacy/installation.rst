Installation
============

To install RascalC, simply clone the Github repository and compile the C++ code (see :ref:`dependencies` below). This is done as follows::

    git clone https://github.com/oliverphilcox/RascalC.git
    cd RascalC
    make

**NB**: RascalC can be run in various modes by adding compiler flags in the ``Makefile``. See :doc:`getting-started` and :doc:`main-code` for more information.

Once RascalC is installed, see the :doc:`getting-started` and :doc:`tutorial` sections.

.. _dependencies:

Dependencies
-------------

RascalC requires the following packages:

- `C compiler <https://gcc.gnu.org/>`_: Tested with 5.4 or later
- `Gnu Scientific Library (GSL) <https://www.gnu.org/software/gsl/doc/html/index.html>`_: Any recent version
- `Corrfunc <https://corrfunc.readthedocs.io>`_: 2.0 or later
- (*Optional but encouraged*) `OpenMP  <https://www.openmp.org/'>`_: Any recent version (required for parallelization)

Corrfunc can be installed using ``pip install corrfunc`` and is used for efficient pair counting.

For the Python pre- and post-processing we require:

- `Python <https://www.python.org/>`_: 2.7 or later, 3.4 or later
- `Numpy <http://www.numpy.org/>`_: 1.10 or later
- (*Optional*) `Healpy <https://healpy.readthedocs.io/en/latest/>`_: any recent version. (Necessary if using HealPix jackknife regions)

These can be simply installed with pip or conda.

.. _acknowledgements:

Acknowledgements
-----------------

Main Authors:

- Oliver H. E. Philcox (Princeton / Harvard)
- Daniel J. Eisenstein (Harvard)
- Ross O'Connell (Pittsburgh)
- Alexander Wiegand (Garching)
- Misha Rashkovetskyi (Harvard)

Please cite `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_, `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_ , `Philcox et al. 2019 <https://arxiv.org/abs/1904.11070>`_ and `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_ when using this code in your research.

Note that many of the code modules and convenience functions are shared with the small-scale power spectrum estimator code `HIPSTER <https://HIPSTER.readthedocs.io>`_, developed by Oliver Philcox and Daniel Eisenstein.
