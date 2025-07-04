Covariance Matrix Estimation
=============================

**NB**: This documentation page is about the legacy/command-line usage. We believe it is much easier to learn the Python library, please see its :doc:`../library/essential-functions` and :doc:`../library/installation`.

Overview
----------

This is the main section of RascalC, where 2PCF or 3PCF covariance matrix estimates are computed via Monte Carlo integration from a given set of input particles. For the 2PCF, depending on the number of input fields the code will compute either components for a single covariance matrix or all required components for 6 cross-covariance matrices (i.e. for multi-tracer analysis).

**Prerequisites**:

- In JACKKNIFE mode (including LEGENDRE_MIX mode with JACKKNIFE), the jackknife weights and binned pair counts must be computed via the :doc:`jackknife-weights` script before the C++ code is run.
- In DEFAULT mode, or LEGENDRE_MIX mode without JACKKNIFE, the RR pair counts must be computed via the :ref:`RR_counts` script before the C++ code is run.
- In LEGENDRE and 3PCF modes, the survey correction functions :math:`\Phi` must be computed via the :ref:`survey_correction_2PCF` script before the C++ code is run.
- In LEGENDRE_MIX mode, the projection coefficients from :math:`\mu` bins to Legendre multipoles musst be computed via the :ref:`mu_bin_legendre_factors` script before the C++ code is run.

.. _particle-grid:

Particle Grid and Cells
~~~~~~~~~~~~~~~~~~~~~~~~~

In the code, the particles are divided up into many cells for efficient computation, with each containing :math:`\sim10` random particles. A cuboidal grid of cells is constructed aligned with the :math:`(x,y,z)` axes, which fully encompasses all random particles. The cellsize is set by the the :math:`N_\mathrm{side}` parameter, which gives the number of (cubic) cells along the largest dimension along the box. This must be an odd integer, and the code will automatically exit if :math:`N_\mathrm{side}` is too small or too large (since this gives inefficient performance).

.. _covariance-precision:

Covariance Matrix Precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The precision of covariance matrix estimators can be user-controlled in the RascalC code. To understand these parameters we must briefly outline the selection of random quads of particles (denoted :math:`\{i,j,k,l\}`) used in the Monte Carlo integration. The particle selection algorithm has the following form:

1. Loop over all :math:`i`-particles :math:`N_\mathrm{loops}` times. Each loop is independent and can be run on separate cores.
2. For each :math:`i` particle, we pick :math:`N_2` :math:`j`-particles at random, according to some selection rule. Here we compute the 2-point contribution to the covariance matrix, if we are computing a 2PCF covariance.
3. For each :math:`j` particle, we pick :math:`N_3` :math:`k`-particles at random, according to some selection rule. Here we compute the 3-point contribution to the covariance matrix.
4. For each :math:`k` particle, we pick :math:`N_4` :math:`l`-particles at random, according to some selection rule. Here we compute the 4-point contribution to the covariance matrix.
5. *(3PCF mode only)*: For each :math:`l` particle, we pick :math:`N_5` :math:`m`-particles at random, according to some selection rule. Here we compute the 5-point contribution to the 3PCF covariance matrix.
6. *(3PCF mode only)*: For each :math:`m` particle, we pick :math:`N_6` :math:`n`-particles at random, according to some selection rule. Here we compute the 6-point contribution to the 3PCF covariance matrix.

By setting the parameters :math:`(N_\mathrm{loops},N_2, N_3, N_4,[N_5,N_6])` we can control the precision of each matrix component. Standard values of :math:`N_2\sim N_3\sim N_4 [\sim N_5 \sim N_6] \sim 10` normally work well. Each loop of the code produces an independent estimate of the full covariance matrix, which can be used to create accurate inverse matrices and effective number of mock calculations. The covariance converges relatively fast, so setting :math:`N_\mathrm{loops}`
to a few times the number of cores should work well. Values of :math:`N_\mathrm{loops}\gtrsim 100` should be avoided to stop file sizes and reconstruction times becoming large.

Note that we require a relatively large value of :math:`N_2N_3N_4N_\mathrm{loops}` for the output matrices to converge, in particular if we wish to invert the matrices. If the output matrices are not sufficiently converged, the reconstruction scripts (:doc:`post-processing`) will fail.

Usage
------

The code is used as follows, with the command line options specified below:

.. code-block:: bash

    make
    ./cov [OPTIONS]

The first line produces the ``./cov`` file by compiling the code; running ``make clean`` is not necessary in most cases since recompilation is invoked automatically when needed depending on changes in source files. The Makefile may need to be altered depending on the particular computational configuration used. The default Makefile is for a standard Unix installation, with the ``Makefile_mac`` file giving a sample Makefile for a Mac installation. This uses the following optional flags in the Makefile;

- ``-DOPENMP``: (Recommended) Run code in parallel with OpenMP, using the OpenMP installation specfied by the ``-lgomp`` and ``-fopenmp`` flags.
- ``-DPERIODIC``: Use periodic boundary conditions (appropriate for a cubic simulation box, but not mock surveys).
- ``-DJACKKNIFE``: Compute both full-survey and jackknife 2PCF covariance matrix terms, allowing for shot-noise-rescaling calibration from the survey itself.
- ``-DLEGENDRE``: Compute the full-survey covariance matrix terms for (even) Legendre multipoles of the 2PCF accumulated directly. Incompatible with jackknives.
- ``-DLEGENDRE_MIX``: Compute the full-survey covariance matrix terms for (even) Legendre multipoles of the 2PCF projected from (a typically large number of) :math:`mu` bins (estimated in this way in `pycorr <https://py2pcf.readthedocs.io>`, for example). Compatible with jackknives; all the counts should be computed with sufficiently large number of :math:`mu` bins, ideally as many of them as the Legendre multipoles are projected from. For jackknife covariance, the disconnected term is dropped (should not make a significant difference since the term has been found tiny in practice).
- ``-DTHREE_PCF``: Compute the full-survey covariance matrix terms for (even and odd) Legendre multipoles of the isotropic 3PCF.
- DEFAULT mode refers to the case when neither ``LEGENDRE`` (nor ``LEGENDRE_MIX``) nor ``THREE_PCF`` are enabled. Then the covariance is computed for :math:`(r,\mu)`-binned correlation function.

**NB**: For a summary of input command line parameters, simply run ``./cov`` with no arguments.

Options
~~~~~~~

Input parameters for the RascalC code may be specified by passing options on the command line or by setting variables in the Parameters class in the ``modules/parameters.h``. When using two sets of random particles, the code will exit automatically if all the required files are not input. We list the major code options below, in the form ``-option`` (*variable*) for command line option ``option`` and Parameters class variable *variable*.

**Essential Parameters**:

- ``-def``: Run the code with the default options for all parameters (as specified in the ``modules/parameters.h`` file.
- ``-in`` (*fname*): Input ASCII random particle file for the first set of tracer particles. This must be in {x,y,z,w,j} format, as described in :ref:`file-inputs`.
- ``-binfile`` (*radial_bin_file*): Radial binning ASCII file (see :ref:`file-inputs`) specifying upper and lower bounds of each radial bin.
- ``-cor`` (*corname*): Input correlation function estimate for the first set of particles in ASCII format, as specified in :ref:`file-inputs`. This can be user defined or created by :ref:`full-correlations`.
- ``-binfile_cf`` (*radial_bin_file_cf*): Radial binning ASCII file for the correlation function (see :ref:`file-inputs`) specifying upper and lower bounds of each radial bin.
- ``-norm`` (*nofznorm*): Number of galaxies in the first set of tracer particles. This is used to rescale the random particle covariances.
- ``-output`` (*out_file*): Output directory in which to store covariance matrix estimates. This directory will be created if not already present. **Beware**: the code can produce a large volume of output (:math:`\sim 1` GB for a standard run with one field and :math:`\sim1000` bins).
- ``-mbin_cf`` (*mbin_cf*): Number of :math:`\mu` bins used for the correlation function.
- ``-nside`` (*nside*): Number of cubic cells to use along the longest dimension of the grid encompassing the random particles, i.e. :math:`N_\mathrm{side}`. See :ref:`particle-grid` note for usage.
- ``-nthread`` (*nthread*): Number of parallel processing threads used if code is compiled with OpenMPI.

**DEFAULT and LEGENDRE_MIX mode Binning Parameters**:

- ``-mbin`` (*mbin*): Number of :math:`\mu` bins used. This must match that used to create the jackknife weights.
- ``-RRbin`` (*RR_bin_file*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}_11.dat`` ASCII file containing the summed pair counts in each bin (:math:`RR_{aA}^{11}`), created by the :file:`jackknife_weights` scripts.

**JACKKNIFE mode Parameters**:

- ``-jackknife`` (*jk_weight_file*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}_11.dat`` file containing the jackknife weights for each bin (:math:`w_{aA}^{11}`), as created by the :file:`jackknife_weights` scripts.

**LEGENDRE and 3PCF mode Parameters**:

- ``max_l`` (*max_l*): Maximum Legendre moment to compute. This must be even in the LEGENDRE or LEGENDRE_MIX mode.
- ``phi_file`` (*phi_file*): Location of the file containing the survey correction function parameters, as created by the :ref:`survey_correction_2PCF` or :ref:`survey_correction_3PCF` script. Must not be given in LEGENDRE_MIX mode.
- ``mu_bin_legendre_file`` (*mu_bin_legendre_file*): Location of the file containing the projection factors from :math:`\mu` bins to Legendre multipoles as produced by the :ref:`mu_bin_legendre_factors` script, only for LEGENDRE_MIX mode.

**Precision Parameters**

- ``-maxloops`` (*max_loops*): This is the number of matrix subsamples to compute. See :ref:`covariance-precision` note for usage guidelines. (Default: 10)
- ``-N2``, ``-N3``, ``-N4`` (*N2*, *N3*, *N4*): The parameters controlling how many random particles to select at each stage. See :ref:`covariance-precision` note above. (Default: 10)
- ``-N5``, ``-N5`` (*N5*, *N6*): As above, but for the 3PCF mode only. (Default: 10)

**General Multi Field Parameters**:

- ``-in2`` (*fname2*): Input ASCII random particle file for the second set of tracer particles.
- (*nofznorm2*): Total number of galaxies in the second set of tracer particles.
- ``-cor12`` (*corname12*): Input cross correlation function file between the two sets of random particles, as created by :ref:`full-correlations`.
- ``-cor2`` (*corname2*): Input autocorrelation function for the second set of particles, either user-defined or created by :ref:`full-correlations`.
- ``-norm2`` (*nofznorm2*): Number of galaxies in the second set of tracer particles. This is used to rescale the random particle covariances.

**DEFAULT and LEGENDRE_MIX mode Multi Field Parameters**:

- ``-jackknife12`` (*jk_weight_file12*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}_12.dat`` file containing the jackknife weights for each bin for the combination of random particle sets 1 and 2 (:math:`w_{aA}^{12}`), as created by the :file:`jackknife_weights` scripts.
- ``-jackknife2`` (*jk_weight_file2*): Location of the ``jackknife_weights_n{N}_m{M}_j{J}_22.dat`` file containing the jackknife weights for each bin for the second set of random particles (:math:`w_{aA}^{22}`), as created by the :file:`jackknife_weights` scripts.
- ``-RRbin12`` (*RR_bin_file12*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}_12.dat`` ASCII file containing the summed jackknife pair counts in each bin for the combination of random particle sets 1 and 2 (:math:`RR_{aA}^{12}`), created by the :file:`jackknife_weights` scripts.
- ``-RRbin2`` (*RR_bin_file2*): Location of the ``binned_pair_counts_n{N}_m{M}_j{J}_22.dat`` ASCII file containing the summed jackknife pair counts in each bin for the combination of random particle sets 1 and 2 (:math:`RR_{aA}^{22}`), created by the :file:`jackknife_weights` scripts.

**LEGENDRE mode Multi Field Parameters**:

- ``phi_file12`` (*phi_file12*): Location of the file containing the survey correction function parameters for the for the second field, as created by the :doc:`geometry-correction` script.
- ``phi_file2`` (*phi_file2*): Location of the file containing the survey correction function parameters for the for the combination of fields 1 and 2, as created by the :doc:`geometry-correction` script.

**Optional Parameters**

- ``-perbox``: This flag signifies a periodic box computation (which requires the code to be compiled with the ``-DPERIODIC`` flag in the Makefile) as opposed to the default aperiodic behavior (which requires the code to be compiled without the ``-DPERIODIC`` flag in the Makefile).
- ``-boxsize`` (*boxsize*): Sets the size of the computational domain if the periodic box is enabled, or the random particles are created in RascalC. (Default: 2000)
- ``-seed`` (*seed*): Random number generator seed. If given, allows to reproduce the results with the same settings, except the number of threads.
- ``-mumin`` (*mumin*): Minimum :math:`\mu` binning to use in the analysis. (Default: 0, or -1 in 3PCF mode)
- ``-mumax`` (*mumax*): Maximum :math:`\mu` binning to use in the analysis. (Default: 1)
- ``-cf_loops`` (*cf_loops*): Number of iterations over which to refine the correlation functions.
- ``-rescale`` (*rescale*): Factor by which to multiply all the input positions. Zero or negative values are reset to the *boxsize*. (Default: 1)
- ``-xicut`` (*xicutoff*): The radius beyond which the correlation functions :math:`\xi(r,\mu)` are set to zero. (Default: 250)
- ``-nmax`` (*nmax*): The maximum number of particles to read in from the random particle files. (Default: 1e12)
- ``-save`` (*savename*): If *savename* is set, the cell selection probability grid is stored as *savename*. This must end in ``.bin``. (Default: NULL)
- ``-load`` (*loadname*): If set, load a cell selection probability grid computed in a previous run of RascalC. (Default: NULL)
- ``-invert`` (*qinvert*): If this flag is passed to RascalC, all input particle weights are multiplied by -1. (Default: 0)
- ``-balance`` (*qbalance*): If this flag is passed to RascalC, all negative weights are rescaled such that the total particle weight is 0. (Default: 0)
- ``-np`` (*np*): If provided, this overrides any input random particle file and creates *np* randomly drawn particles in the cubic box with the side length *boxsize*. **NB**: Currently creating particles at random is only supported for a single set of tracer particles.
- ``-rs`` (*rstart*): If inverting particle weights, this sets the index from which to start weight inversion. (Default: 0)

.. _code-output:

Code Output
-----------

In the specified output directory, RascalC creates the directories ``3PCFCovMatricesAll/`` (3PCF mode), ``CovMatricesAll/`` (DEFAULT, LEGENDRE and JACKKNIFE modes) and ``CovMatricesJack/`` (JACKKNIFE mode) containing the relevant output matrix estimates. These contain multiple estimates of the each part of the total matrix and should be reconstructed using the :doc:`post-processing` scripts.

The full output files take the following form (for N radial bins, M angular bins, maximum Legendre bin L and J non-zero jackknife regions, with FIELDS specifying the utilized tracer fields):

*3PCF or LEGENDRE mode*:

- ``c{X}_n{N}_l{L}_{FIELDS}_{I}.txt``: I-th estimate of the X-point covariance matrix estimates, i.e. :math:`{}^X\mathbf{C}`. The summed covariance matrix terms have the suffix 'full'.
- ``binct_c{X}_n{N}_l{L}_{FIELDS}_{I}.txt``: Total used counts per bin for the X-point covariance matrix.
- ``total_counts_n{N}_l{L}_{FIELDS}_{I}``: Total number of sets of particles attempted for the summed integral.

*DEFAULT or JACKKNIFE mode*:

 - ``c{X}_n{N}_m{M}_j{J}_{FIELDS}_{I}.txt``: I-th estimate of the X-point covariance matrix estimates, i.e. :math:`C_{X,ab}` The summed covariance matrix has the suffix 'full'.
 - ``RR_n{N}_m{M}_{FIELDS}_{I}.txt``: I-th estimate of the (non-jackknife) :math:`RR_{ab}^{XY}` pair counts which can be compared with Corrfunc.
 - ``binct_c{X}_n{N}_m{M}_{FIELDS}.txt``: Total used counts per bin for the X-point covariance matrix.
 - ``total_counts_n{N}_m{M}_{FIELDS}.txt``: Total number of pairs, triples and quads attempted for the summed integral.

 *JACKKNIFE mode only*:

 - ``RR{P}_n{N}_m{M}_{FIELDS}.txt``: Estimate of :math:`RR_{ab}` pair count for particles in random-subset P (:math:`P\in[1,2]`).  This is used to compute the disconnected jackknife matrix term.
 - ``EE{P}_n{N}_m{M}_{FIELDS}.txt``: Estimate of :math:`EE_{ab}` :math:`\xi`-weighted pair count for particles in random-subset P. This is also used for the disconnected jackknife matrix term.

Each file is an ASCII format file containing the relevant matrices with the collapsed bin indices :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` (2PCF) or :math:`\mathrm{bin}_\mathrm{collapsed} = \left(\mathrm{bin}_\mathrm{radial,1}\times n_r + \mathrm{bin}_\mathrm{radial,2}\right)\times n_\mu + \mathrm{bin}_\mathrm{angular}` (3PCF) for a total of :math:`n_\mu` angular (or Legendre) bins and :math:`n_r` radial bins.
