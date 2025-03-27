Package installation
====================

.. highlight:: bash

Simplified way for DESI members at NERSC
----------------------------------------

Recommended to use with `cosmodesi` environment.
In particular, load it before installing::

    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
    pip install -e /global/common/software/desi/users/mrash/RascalC

This installs the library from a common software folder in the development mode, so that after I update it e.g. with some fix, you will have the new version without the need to re-install.

Generic installation
--------------------

Prerequisites
^^^^^^^^^^^^^

OS: Linux (recommended) or macOS. Others (including Windows) not supported.

-  C and C++ compilers: `GNU <https://gcc.gnu.org/>`_ (Linux) or Apple ``clang`` (macOS);
- `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/doc/html/index.html>`_;
- `OpenMP  <https://www.openmp.org/'>`_ (required for parallelization).

On Linux, install these from your distribution's repositories.

On macOS, install the following with `Homebrew <https://brew.sh/>`_ (other setups will likely require editing the ``Makefile``)::

    brew install gsl pkg-config libomp

Additionally, the code requires `pycorr <https://github.com/cosmodesi/pycorr>`_ to deal with pair counts and data correlation function estimators.
To compute pair counts of catalogs, you need a `custom version of Corrfunc <https://github.com/adematti/Corrfunc>`_ (see also `pycorr installation instructions <https://py2pcf.readthedocs.io/en/latest/user/building.html>`_).
Both can be installed with a single command::

    python3 -m pip install 'git+https://github.com/cosmodesi/pycorr#egg=pycorr[corrfunc]'

If you get the GPU support problem with Corrfunc (``Error: To compile with GPU support define "CUDA_HOME" Else set "USE_GPU=0"``), it is probably easiest to prepend either definition to the command: ``CUDA_HOME=... python3 -m pip install 'git+https://github.com/cosmodesi/pycorr#egg=pycorr[corrfunc]'`` or ``USE_GPU=0 python3 -m pip install 'git+https://github.com/cosmodesi/pycorr#egg=pycorr[corrfunc]'``.
Alternatively, you can ``export CUDA_HOME=...`` or ``export USE_GPU=0`` in your shell before running ``python3 -m pip install 'git+https://github.com/cosmodesi/pycorr#egg=pycorr[corrfunc]'``.

One of the reasons we recommend Linux is that building Corrfunc with multi-threading support on macOS has been a very hard experience.

RascalC Python package
^^^^^^^^^^^^^^^^^^^^^^

Clone the Github repository and install the package from the source directory::

    git clone https://github.com/oliverphilcox/RascalC
    cd RascalC
    python3 -m pip install .

Make sure to reinstall after pulling updates.

Or you can install with a single command::

    python3 -m pip install https://github.com/oliverphilcox/RascalC.git

Repeat this command to update the package.