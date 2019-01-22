
Correlation Functions
=======================

The scripts described below are wrappers of the `Corrfunc <https://corrfunc.readthedocs.io>`_ code (Sinha & Garrison 2017), used to create full-survey and jackknife correlation functions. The former are used in the computation of the Gaussian covariance matrices, and the latter allow for determination of the shot-noise rescaling parameter. If the correlation function is required to be computed in a different manner, user-input correlation functions can simply replace the output of these codes, with the file-types described in :ref:`file-inputs`.


.. _full-correlations:

Full Matrix Correlations :math:`\xi(r,\mu)`
----------------------------------------------

To compute the covariance matrix estimates :math:`\hat{C}_{ab}` we require some estimate of the correlation function. Here, we compute the full-survey correlation function with a specified binning using Corrfunc. We provide routines for both 1- and 2-field scenarios (computing the fields :math:`\{\xi^{XX}(r,\mu), \xi^{XY}(r,\mu), \xi^{YY}(r,\mu)\}` in the latter case). This uses both the galaxies and random particle files, and requires a correlation function binning file, such as created by :ref:`write-binning-file`. The estimations are computed via the Landy-Szalay estimator using :math:`\xi^{XY}_a = (DD_a^{XY} - DR_a^{XY} - DR_a^{YX} + RR_a^{XY})/RR_a^{XY}` for bin :math:`a`, fields :math:`X, Y` with DD/DR/RR specifying data-data / data-random / random-random pair counts. 

.. todo:: add note on xi interpolation / refining

*Periodicity*: This script can be run for periodic or aperiodic input data; this corresponds to measuring :math:`\mu` from the :math:`z` or radial axis respectively. If the data is periodic (e.g. from a cosmological box simulation) the -DPERIODIC flag should be set on compilation of the full C++ code in :doc:`main-code`.

**Usage**

For a single field analysis::

    python python/xi_estimator.py {GALAXY_FILE} {RANDOM_FILE} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts}]
    

For an analysis using two distinct fields::
    
    python python/xi_estimator_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1} {RANDOM_FILE_2} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts_11} {RR_counts_12} {RR_counts_22}]

**NB**: The two field script computes all three distinct (cross-)correlations between the two random fields, thus has a runtime :math:`\sim` 3 times that of ``xi_estimator.py``. The two fields should be distinct to avoid issues with double counting. 


**Input Parameters**

- {GALAXY_FILE}, {GALAXY_FILE_1}, {GALAXY_FILE_2}: Input ASCII file containing galaxy positions and weights in {x,y,z,weight,jackknife_ID} format such as that created with the :doc:`pre-processing` scripts.  (Jackknives are not used in this script and may be omitted). This should be in ``.csv``, ``.txt`` or ``.dat`` format with space-separated columns.
- {RANDOM_FILE}, {RANDOM_FILE_1}, {RANDOM_FILE_2}: Input ASCII file containing random particle positions and weights, as for the galaxy files.
- {RADIAL_BIN_FILE}: ASCII file specifying the radial bins for :math:`\xi(r,\mu)`, as described in :ref:`file-inputs`. This can be user-defined or created by the :ref:`write-binning-file` scripts.  **NB**: This bin-file specifies the bins for the *correlation function*, which may be distinct from the *covariance-matrix* bins. In particular, the lowest bin should extend to :math:`r = 0`.
- {MU_MAX}: Maximum :math:`\mu = \cos\theta` used in the angular binning.
- {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu]`.
- {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
- {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
- {OUTPUT_DIR}: Directory in which to house the correlation functions. This will be created if not in existence.
- *Optional* {RR_counts}, {RR_counts_XY}: Pre-computed RR pair counts (for the single field and between fields X and Y respectively). These should be in the format described in :ref:`file-inputs`, and must use the same number of radial and angular bins as specified above. If not specified, these are recomputed by the code. 


**Output Files**

ASCII files are created specifying the correlation function in the file-format given in :ref:`file-inputs`. The filename has the format ``xi_n{N}_m{M}_{INDEX}.dat``, where N and M specify the number of radial and angular bins respectively. INDEX specifies the correlation function type, where 11 = field 1 auto-correlation, 22 = field 2 auto-correlation, 12 = cross-correlation of fields 1 and 2. The first and second lines of the ``.dat`` file list the radial and angular bin centers, then each subsequent line lists the :math:`\xi(r,\mu)` estimate, with the column specifying the :math:`\mu` bin and the row specifying the :math:`r` bin.

    
.. _jackknife-correlations:

Jackknife Matrix Correlations :math:`\xi^J(r,\mu)`
----------------------------------------------------

For later comparison of the jackknife covariance matrix estimate with the data, we require the jackknife covariance matrix, which is derived from the correlation function estimates in each unrestricted jackknife. The scripts below are provided to compute these using Corrfunc. For jackknife :math:`J` and fields :math:`\{X,Y\}`, we compute the pair counts :math:`FG^{XY}_a` in bin :math:`a` (where :math:`F,G\in[D,R]` for data and random fields D and R), from a cross-pair counts between particles in jackknife :math:`A` of :math:`F^X` and the entire of field :math:`G^Y`. These are added to the pair counts from the cross of particles in jackknife :math:`A` of field :math:`G^Y` with the entire of field :math:`F^X` if the fields are distinct. This allows us to compute all :math:`n_\mathrm{jack}` correlation functions :math:`\xi^{XY}_A(r,\mu)` via the Landy-Szalay estimator :math:`\xi^{XY}_{aA} = (DD_{aA}^{XY} - DR_{aA}^{XY} - DR_{aA}^{YX} + RR_{aA}^{XY})/RR_{aA}^{XY}` for bin :math:`a`

**NB**: The binning file used here should be the same as that used for the *covariance matrix* **not** the full correlation function, to allow comparison with the :math:`C^J_{ab}` estimate.


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
