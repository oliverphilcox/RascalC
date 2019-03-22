Package Installation
=====================

To install RascalC, simply clone the Github repository and compile the C++ code (see :ref:`dependencies` below). This is done as follows::

    git clone https://github.com/oliverphilcox/RascalC.git
    cd RascalC
    make

Once RascalC is installed, see the :doc:`getting-started` and :doc:`tutorial` sections.

.. _dependencies:

Dependencies
-------------

RascalC requires the following packages:

- `C compiler <https://gcc.gnu.org/>`_: Tested with gcc 5.4.O
- `Gnu Scientific Library (GSL) <https://www.gnu.org/software/gsl/doc/html/index.html>`_: Any recent version
- `Corrfunc <https://corrfunc.readthedocs.io>`_: 2.0 or later
- (*Optional but encouraged*) `OpenMP  <https://www.openmp.org/'>`_: Any recent version (required for parallelization)

Corrfunc can be installed using ``pip install corrfunc`` and is used for efficient pair counting.

For the Python pre- and post-processing we require:

- `Python <https://www.python.org/>`_: 2.7 or later, 3.4 or later
- `Numpy <http://www.numpy.org/>`_: 1.10 or later
- (*Optional*) `Healpy <https://healpy.readthedocs.io/en/latest/>`_: any recent version. (Necessary if using HealPix jackknife regions)

These can be simply installed with pip or conda.

