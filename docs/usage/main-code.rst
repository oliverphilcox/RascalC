Covariance Matrix Estimation
=============================

Overview
----------

This is the main section of RascalC, where a covariance matrix estimates are computed via Monte Carlo integration from a given set of input particles. Depending on the number of input fields the code will compute either components for a single covariance matrix or all required components for 6 cross-covariance matrices. 

**NB**: Before running this code, the jackknife weights and binned pair counts must be computed 

.. _particle-grid:

Particle Grid and Cells
~~~~~~~~~~~~~~~~~~~~~~~~~

In the code, the particles are divided up into many cells for efficient computation, with each containing :math:`\sim10` random particles. A cuboidal grid of cells is constructed aligned with the :math:`(x,y,z)` axes, which fully encompasses all random particles. The cellsize is set by the the *nside* parameter, which gives the number of (cubic) cells along the largest dimension along the box. This must be an odd integer, and the code will automatically exit if *nside* is too small or too large (since this gives inefficient performance).

.. _covariance-precision:

Covariance Matrix Precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The precision of covariance matrix estimators can be user-controlled in the RascalC code. To understand these parameters we must briefly outline the selection of random quads of particles (denoted :math:`\{i,j,k,l\}`) used in the Monte Carlo integration. The particle selection algorithm has the following form:

1. Loop over all :math:`i` particles :math:`N_\mathrm{loops}` times. Each loop is independent and can be run on separate cores.
2. For each :math:`i` particle, we pick :math:`N_2` :math:`j` particles at random, according to some selection rule. Here we compute the 2-point contribution to the covariance matrix.
3. For each :math:`j` particle, we pick :math:`N_3` :math:`k` particles at random, according to some selection rule. Here we compute the 3-point contribution to the covariance matrix.
4. For each :math:`k` particle, we pick :math:`N_4` :math:`l` particles at random, according to some selection rule. Here we compute the 4-point contribution to the covariance matrix.

By setting the parameters :math:`(N_\mathrm{loops},N_2, N_3, N_4)` we can control the precision of each matrix component. Standard values of :math:`N_2\sim N_3\sim N_4 \sim 10` normally work well. Each loop of the code produces an independent estimate of the full covariance matrix, which can be used to create accurate inverse matrices and effective number of mock calculations. The covariance converges relatively fast, so setting :math:`N_\mathrm{loops}` 
to a few times the number of cores should work well.

Usage
------

The code is used as follows, with the command line options specified below::
    
    bash clean
    make 
    ./cov [OPTIONS]

The first line removes any pre-existing C++ file before it is recompiled in line 2 to produce the ``./cov`` file. For a summary of input command line parameters, simply run ``./cov`` with no arguments.

**NB**: For periodic random particle files (such as those produced by simulations), the code should be compiled with the ``-DPERIODIC`` flag in the Makefile. This is turned off by default.

Options
~~~~~~~

Input parameters for the RascalC code may be specified by passing options on the command line or by setting variables in the Parameters class in the ``modules/parameters.h``. When using two sets of random particles, the code will exit automatically if all the required files are not input. We list the major code options below, in the form ``-option`` (*variable*) for command line option ``option`` and Parameters class variable *variable*.

**Essential Parameters**:

- ``-def``, ``-default``: Run the code with the default options for all parameters (as specified in the ``modules/parameters.h`` file.
- ``-in``, ``-in1`` (*fname*): Input ASCII random particle file for the first set of tracer particles. This must be in {x,y,z,w,j} format, as described in :ref:`file-inputs`.
- ``-binfile``, ``-radial_binfile`` (*radial_bin_file*): Radial binning ASCII file (see :ref:`file-inputs`) specifying upper and lower bounds of each radial bin.
- ``-cor``, ``-cor1`` (*corname*): Input correlation function estimate for the first set of particles in ASCII format, as specified in :ref:`file-inputs`.
- ``-norm`` (*nofznorm*): Total number of galaxies in the first set of tracer particles. This is used to rescale the random particle covariances.
- ``-jackknife``, ``-jackknife_weights``, ``-jackknife1``, ``-jackknife_weights1`` (*jk_weight_file*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}.dat`` file containing the jackknife weights for each bin (:math:`w_{aA}^{11}`), as created by the :file:`jackknife_weights` scripts.
- ``-RRbin``, ``-RRbin_file``, ``-RRbin1``, ``RRbin_file1`` (*RR_bin_file*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}.dat`` ASCII file containing the summed jackknife pair counts in each bin (:math:`RR_{aA}^{11}`), created by the :file:`jackknife_weights` scripts.
- ``-output`` (*out_file*): Output directory in which to store covariance matrix estimates. This directory will be created if not already present. **Beware**: the code can produce a large volume of output (:math:`\sim 1` GB for a standard run with one field and :math:`\sim1000` bins). 
- ``-mbin`` (*mbin*): Number of :math:`\mu` bins used. This must match that used to create the jackknife weights. 
- ``-nthread``, ``-nthreads`` (*nthread*): Number of parallel processing threads used if code is compiled with OpenMPI.
- ``-nside``, ``-ngrid``, ``-grid`` (*nside*): Number of cubic cells to use along the longest dimension of the grid encompassing the random particles. See :ref:`particle-grid` note for usage.

**Additional Multi Field Parameters**:

- ``-in2`` (*fname2*): Input ASCII random particle file for the second set of tracer particles.
- (*nofznorm2*): Total number of galaxies in the second set of tracer particles.
- ``-cor12`` (*corname12*): Input cross correlation function file between the two sets of random particles, as created by **XXX**.
- ``-cor2`` (*corname2*): Input autocorrelation function for the second set of particles, either user-defined or created by **XXX**.

.. todo:: add in correlation function creator script

- ``-jackknife12``, ``-jackknife_weights12`` (*jk_weight_file12*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}.dat`` file containing the jackknife weights for each bin for the combination of random particle sets 1 and 2 (:math:`w_{aA}^{12}`), as created by the :file:`jackknife_weights` scripts.
- ``-jackknife2``, ``-jackknife_weights2`` (*jk_weight_file12*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}.dat`` file containing the jackknife weights for each bin for the second set of random particles (:math:`w_{aA}^{22}`), as created by the :file:`jackknife_weights` scripts.
- ``-RRbin12``, ``-RRbin_file12`` (*RR_bin_file12*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}.dat`` ASCII file containing the summed jackknife pair counts in each bin for the combination of random particle sets 1 and 2 (:math:`RR_{aA}^{12}`), created by the :file:`jackknife_weights` scripts.
- ``-RRbin2``, ``-RRbin_file2`` (*RR_bin_file2*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}.dat`` ASCII file containing the summed jackknife pair counts in each bin for the combination of random particle sets 1 and 2 (:math:`RR_{aA}^{22}`), created by the :file:`jackknife_weights` scripts.

**Precision Parameters**

- ``-maxloops``, ``-loops`` (*max_loops*): This is the number HEDWSDHEJIL See :ref:`covariance-precision` note for usage guidelines. (Default: 10)
- (*N2*, *N3*, *N4*): The parameters controlling how many random particles to select at each stage. See :ref:`covariance-precision` note above. (Default: 10)

**Optional Parameters**

- ``-mumin`` (*mumin*): Minimum :math:`\mu` binning to use in the analysis. (Default: 0) 
- ``-mumax`` (*mumax*): Maximum :math:`\mu` binning to use in the analysis. (Default: 1)
- (*perbox*): Boolean controlling whether we are using a periodic box. (Default: False)
- ``-boxsize`` (*boxsize*): If creating particles randomly, this is the periodic size of the computational domain. If particles are read from file, this is set dynamically. (Default: 400)
- ``-rescale``, ``-scale`` (*rescale*): Factor by which to dilate the input positions. Zero or negative values cause this to be set to the boxsize. (Default: 1)
- ``-xicut`` (*xicutoff*): The radius beyond which the correlation functions :math:`\xi(r,\mu)` are set to zero. (Default: 400)
- ``-rs`` (*rstart*):
- ``-nmax`` (*nmax*):
- ``-save``, ``-store`` (*savename*): If *savename* is set, the cell selection probability grid is stored as *savename*. This must end in ``.bin``. (Default: NULL)
- ``-load`` (*loadname*): If set, load a cell selection probability grid computed in a previous run of RascalC. (Default: NULL) 
- ``-invert`` (*qinvert*): If this flag is passed to RascalC, all input particle weights are multiplied by -1. (Default: 0)
- ``-balance`` (*qbalance*): If this flag is passed to RascalC, all negative weights are rescaled such that the total particle weight is 0. (Default: 0)
- ``-ran``, ``-np`` (*np*, *make_random*): If *make_random*=1, this overrides any input random particle file and creates *np* randomly drawn particles in the cubic box. **NB**: The command line argument automatically sets *make_random* = 1. 
- (*rstart*): If inverting particle weights, this sets the index from which to start weight inversion. (Default: 0)

.. todo:: don't let code run with random particle creation and multiple fields. And note this in doc somewhere.

.. todo:: add in N2, N3, N4 as free parameters. Remove some of these limits? Add in nofznorm2 as input. Add rstart, perbox as input.


.. _code-output:

Code Output
-----------

In the specified output directory, RascalC creates two directories; ``CovMatricesAll/`` and ``CovMatricesJack`` containing total and jackknife covariance matrix estimates respectively. These contain multiple estimates of the each part of the total matrix and should be reconstructed using the :doc:`post-processing` scripts.

The full output files take the following form (for N radial bins, M radial bins and J non-zero jackknife regions, with FIELDS specifying the utilized tracer fields):
 - ``c{X}_n{N}_m{M}_j{J}_{FIELDS}_{I}.txt``: I-th estimate of the X-point covariance matrix estimates, i.e. :math:`C_{X,ab}` The summed covariance matrix has the suffix 'full'. 
 - ``RR_n{N}_m{M}_{FIELDS}_{I}.txt``: I-th estimate of the (non-jackknife) :math:`RR_{ab}^{XY}` pair counts which can be compared with Corrfunc.
 - ``binct_c{X}_n{N}_m{M}_{FIELDS}.txt``: Total used counts per bin for the X-point covariance matrix.
 - ``total_counts_n{N}_m{M}_{FIELDS}.txt``: Total number of pairs, triples and quads attempted for the summed integral.
 - ``RR{P}_n{N}_m{M}_{FIELDS}.txt``: Estimate of :math:`RR_{ab}` pair count for particles in random-subset P (:math:`P\in[1,2]`).  This is used to compute the disconnected jackknife matrix term.
 - ``EE{P}_n{N}_m{M}_{FIELDS}.txt``: Estimate of :math:`EE_{ab}` :math:`\xi`-weighted pair count for particles in random-subset P. This is also used for the disconnected jackknife matrix term.

Each file is an ASCII format file containing the relevant matrices with the collapsed bin indices :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins. 

