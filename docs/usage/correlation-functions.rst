
Correlation Functions
=======================

The scripts described below are wrappers of the `Corrfunc <https://corrfunc.readthedocs.io>`_ code (Sinha & Garrison 2017), used to create full-survey and jackknife correlation functions. The former are used in the computation of the Gaussian covariance matrices, and the latter allow for determination of the shot-noise rescaling parameter. If the correlation function is required to be computed in a different manner, user-input correlation functions can simply replace the output of these codes, with the file-types described in :ref:`file-inputs`.


.. _full-correlations:

Full Matrix Correlations :math:`\xi(r,\mu)`
----------------------------------------------


**Usage**

For a single field analysis::

    python python/XXXX
    

For an analysis using two distinct fields::
    
    python python/XXXX
    

**Output Files**

    
.. _jackknife-correlations:

Jackknife Matrix Correlations :math:`\xi^J(r,\mu)`
----------------------------------------------------

**Usage**

For a single field analysis::

    python python/XXXX


For an analysis using two distinct fields::
    
    python python/XXXX

    
**Output Files**


.. 
.. 
.. 
.. Here, we compute the weights assigned to each jackknife region for each bin. This is done using the `Corrfunc <https://corrfunc.readthedocs.io>`_ of Sinha & Garrison to compute the weights :math:`w_{aA}^{XY} = RR_{aA}^{XY} / \sum_B RR_{aB}^{XY}` for bin :math:`a`, jackknife :math:`A` and fields :math:`X` and :math:`Y`. 
.. 
.. Two codes are supplied; one using a single set of tracer particles and the other with two input sets, for computation of cross-covariance matrices. These are in the ``python/`` directory. This must be run before the main C++ code.
.. 
.. Usage
.. ~~~~~~~
.. For a single field analysis::
.. 
..     python python/jackknife_weights.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} OUTPUT_DIR}
.. 
.. For an analysis using two distinct fields::
.. 
..     python python/jackknife_weights_cross.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR}
..     
.. **NB**: The two field script computes all three combinations of weights between the two random fields, thus has a runtime :math:`\sim` 3 times that of ``jackknife_weights.py``. Running these together in one script ensures that we have the same number of jackknives for all fields. Also, the two fields must be distinct, else there are issues with double counting. 
.. 
.. .. todo:: check RascalC read-in procedure with all weights 
..     
.. **Input Parameters**
.. 
.. - {RANDOM_PARTICLE_FILE}, {RANDOM_PARTICLE_FILE_1}, {RANDOM_PARTICLE_FILE_2}: Input ASCII file containing random particle positions and jackknife numbers in {x,y,z,weight,jackknife_ID} format, such as that created with the :doc:`pre-processing` scripts. This should be in ``.csv``, ``.txt`` or ``.dat`` format with space-separated columns.
.. - {BIN_FILE}: ASCII file specifying the radial bins, as described in :ref:`file-inputs`. This can be user-defined or created by the :ref:`write-binning-file` scripts.
.. - {MU_MAX}: Maximum :math:`\mu = \cos\theta` used in the angular binning.
.. - {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu]`.
.. - {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
.. - {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
.. - {OUTPUT_DIR}: Directory in which to house the jackknife weights and pair counts. This will be created if not in existence.
.. 
.. 
.. **Notes**:
.. 
.. - This is a very CPU intensive computation since we must compute pair counts between every pair of random particles. The process can be expedited using multiple CPU cores or a reduced number of random particles (e.g. via the :ref:`particle-subset` script).
.. - For two sets of input particles, three sets of weights must be computed for the three possible pairs of two distinct fields, hence the computation time increases by a factor of three.
.. 
.. **Note on Periodicity**
.. 
.. The code can be run for datasets created with either periodic or non-periodic boundary conditions. Periodic boundary conditions are often found in cosmological simlulations. If periodic, the pair-separation angle :math:`\theta` (used in :math:`\mu=\cos\theta`) is measured from the :math:`z` axis, else it is measured from the radial direction. If periodic data is used, the C++ code **must** be compiled with the -DPERIODIC flag.
.. 
.. Output files
.. ~~~~~~~~~~~~~
.. 
.. This code creates ASCII files containing the jackknife weights for each bin and the RR pair counts. The output files have the format ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` where N and M specify the number of radial and angular bins respectively and J gives the number of non-empty jackknives. INDEX specifies which fields are being used i.e. INDEX = 12 implies the :math:`w_{aA}^{12}` and :math:`RR_a^{12}` quantities.
.. 
.. The binned pair counts is a list of weighted pair counts for each bin, summed over all jackknife regions, in the form :math:`RR_a^{J,XY} = \sum_B RR_{aB}^{XY}`. The jackknife weights file lists the weights :math:`w_{aA}^{XY}` for each bin and jackknife region. The :math:`j`-th row contains the (tab-separated) weights for each bin using the :math:`j`-th jackknife. The first value in each row is the jackknife number, and the bins are ordered using the collapsed binning :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins.  
