
Correlation Functions
=======================

The scripts described below are wrappers of the `Corrfunc <https://corrfunc.readthedocs.io>`_ code (Sinha & Garrison 2017), used to create full-survey and jackknife correlation functions. The former are used in the computation of the Gaussian covariance matrices, and the latter allow for determination of the shot-noise rescaling parameter. If the correlation function is required to be computed in a different manner, user-input correlation functions can simply replace the output of these codes, with the file-types described in :ref:`file-inputs`.

.. _full-correlations:

Full Survey Correlation Functions 
-----------------------------------

To compute the covariance matrix estimates :math:`\hat{C}_{ab}` we require some estimate of the correlation function. Here, we compute the full-survey correlation function with a specified binning using Corrfunc. We provide routines for both 1- and 2-field scenarios (computing the fields :math:`\{\xi^{11}(r,\mu), \xi^{12}(r,\mu), \xi^{22}(r,\mu)\}` in the latter case). This uses both the galaxies and random particle files, and requires a correlation function binning file, such as created by :ref:`write-binning-file`. The estimations are computed via the Landy-Szalay estimator using :math:`\xi^{XY}_a = (\widehat{DD}_a^{XY} - \widehat{DR}_a^{XY} - \widehat{DR}_a^{YX} + \widehat{RR}_a^{XY})/\widehat{RR}_a^{XY}` for bin :math:`a`, fields :math:`X, Y` with DD/DR/RR specifying data-data / data-random / random-random pair counts. The hats indicate that the quantities are normalized by the product of the summed weights in the two fields (i.e. :math:`\sum_{i\in X}w_i^X\sum_{j\in Y}w_j^Y`). Note that the binned correlation function estimates :math:`\hat\xi^{XY}_a` cannot simply placed at the bin-centers and interpolated to give a smooth :math:`\xi^{XY}(r,\mu)` estimate; here we use an iterative approach to convert these estimates into interpolation points for the smooth :math:`\xi^{XY}(r,\mu)`.

The scripts take two sets of random particle files; one to compute :math:`DR` counts and one to compute :math:`RR` counts. This allows for a larger number of randoms to be used for :math:`RR` counts, as is often useful.

*Periodicity*: This script can be run for periodic or aperiodic input data; this corresponds to measuring :math:`\mu` from the :math:`z` or radial axis respectively. If the data is periodic (e.g. from a cosmological box simulation) the -DPERIODIC flag should be set on compilation of the full C++ code in :doc:`main-code`.

**Usage**

For a single field analysis::

    python python/xi_estimator.py {GALAXY_FILE} {RANDOM_FILE} {RADIAL_BIN_FILE_DR} {RADIAL_BIN_FILE_RR} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts}]
    

For an analysis using two distinct fields::
    
    python python/xi_estimator_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1_DR} {RANDOM_FILE_1_RR} {RANDOM_FILE_2_DR} {RANDOM_FILE_2_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts_11} {RR_counts_12} {RR_counts_22}]

**NB**: The two field script computes all three distinct (cross-)correlations between the two fields, thus has a runtime :math:`\sim` 3 times that of ``xi_estimator.py``. The two fields should be distinct to avoid issues with double counting. 


**Input Parameters**

- {GALAXY_FILE}, {GALAXY_FILE_1}, {GALAXY_FILE_2}: Input ASCII file containing galaxy positions and weights in {x,y,z,weight,jackknife_ID} format such as that created with the :doc:`pre-processing` scripts.  (Jackknives are not used in this script and may be omitted). This should be in ``.csv``, ``.txt`` or ``.dat`` format with space-separated columns.
- {RANDOM_FILE_DR}, {RANDOM_FILE_1_DR}, {RANDOM_FILE_2_DR}: Input ASCII file containing random particle positions and weights to be used for DR pair counting (with filetype as for the galaxy files).
- {RANDOM_FILE_RR}, {RANDOM_FILE_1_RR}, {RANDOM_FILE_2_RR}: Input ASCII file containing random particle positions and weights to be used for RR pair counting (with filetype as for the galaxy files). **NB**: If pre-computed RR pair counts are specified only the length of the RR random file is used by the code (for normalization).
- {RADIAL_BIN_FILE}: ASCII file specifying the radial bins for :math:`\xi(r,\mu)`, as described in :ref:`file-inputs`. This can be user-defined or created by the :ref:`write-binning-file` scripts.  **NB**: This bin-file specifies the bins for the *correlation function*, which may be distinct from the *covariance-matrix* bins. In particular, the lowest bin should extend to :math:`r = 0`.
- {MU_MAX}: Maximum :math:`\mu = \cos\theta` used in the angular binning.
- {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu_\mathrm{max}]`.
- {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
- {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
- {OUTPUT_DIR}: Directory in which to house the correlation functions. This will be created if not in existence.
- *Optional* {RR_counts}, {RR_counts_XY}: Pre-computed RR pair counts (for the single field and between fields X and Y respectively). These should be in the format described in :ref:`file-inputs`, and must use the same number of radial and angular bins as specified above. If not specified, these are recomputed by the code. Since the full correlation function typically uses a different binning to the output covariance matrix, we typically cannot use the pair counts computed in :doc:`jackknife-weights` and must recompute them. In addition, these should be normalized by the squared sum of weights :math:`(\sum_i w_i)^2` where :math:`i` runs across all random particles in the dataset.


**Output Files**

ASCII files are created specifying the correlation function in the file-format given in :ref:`file-inputs`. The filename has the format ``xi_n{N}_m{M}_{INDEX}.dat``, where N and M specify the number of radial and angular bins respectively. INDEX specifies the correlation function type, where 11 = field 1 auto-correlation, 22 = field 2 auto-correlation, 12 = cross-correlation of fields 1 and 2. The first and second lines of the ``.dat`` file list the radial and angular bin centers, then each subsequent line lists the :math:`\xi(r,\mu)` estimate, with the column specifying the :math:`\mu` bin and the row specifying the :math:`r` bin.

**NB**: The code also prints the number of galaxies in each dataset to the terminal, :math:`N_\mathrm{gal}`. This quantity is important for later normalization of the C++ code.
    
.. _jackknife-correlations:

Jackknife Matrix Correlation Functions
----------------------------------------------

For later comparison of the jackknife covariance matrix estimate with the data, we require the jackknife covariance matrix, which is derived from the correlation function estimates in each unrestricted jackknife. The scripts below are provided to compute these using Corrfunc. For jackknife :math:`J` and fields :math:`\{X,Y\}`, we compute the pair counts :math:`FG^{XY}_a` in bin :math:`a` (where :math:`F,G\in[D,R]` for data and random fields D and R), from a cross-pair counts between particles in jackknife :math:`A` of :math:`F^X` and the entire of field :math:`G^Y`. These are added to the pair counts from the cross of particles in jackknife :math:`A` of field :math:`G^Y` with the entire of field :math:`F^X` if the fields are distinct. This allows us to compute all :math:`n_\mathrm{jack}` correlation functions :math:`\xi^{XY}_A(r,\mu)` via the Landy-Szalay estimator :math:`\xi^{XY}_{aA} = (\widehat{DD}_{aA}^{XY} - \widehat{DR}_{aA}^{XY} - \widehat{DR}_{aA}^{YX} + \widehat{RR}_{aA}^{XY})/\widehat{RR}_{aA}^{XY}` for bin :math:`a`. As before, the code takes two random particle fields of each type, allowing different sized random fields to be used for DR and RR pair counting. For convenience the quantities are normalized by the summed weights across the **entire** set of particles, not just those specific to the given jackknife. The jackknife correlation functions are thus not quite true estimates of :math:`\xi_a`, since they neglect differences in the ratio of galaxies and random particles between galaxies. 

**NB**: The binning file used here should be the same as that used for the *covariance matrix* **not** the full correlation function, to allow comparison with the :math:`C^J_{ab}` estimate.

**Usage**

For a single field analysis::

    python python/xi_estimator_jack.py {GALAXY_FILE} {RANDOM_FILE_DR} {RANDOM_FILE_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_jackknife_counts}]


For an analysis using two distinct fields::
    
    python python/xi_estimator_jack_cross.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1_DR} {RANDOM_FILE_1_RR} {RANDOM_FILE_2_DR} {RANDOM_FILE_2_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_jackknife_counts_11} {RR_jackknife_counts_12} {RR_jackknife_counts_22}]

    
This computes estimates of the auto- and cross-correlations for all unrestricted jackknife regions. Since there are three distinct correlations for each, the run-time is increased by a factor of 3.

Following computation of :math:`\xi^J_{aA}` we can estimate the single-survey jackknife covariance matrix via :math:`C^J_{ab,\mathrm{data}} = \sum_A w_{aA}w_{bA}(\xi^J_{aA}-\bar{\xi}^J_a)(\xi^J_{bA}-\bar{\xi}^J_b) / (1-\sum_B w_{aB}w_{bB})`. This is done internally in the :ref:`post-processing-single` code.

**Input Parameters**

See the input parameters for the :ref:`full-correlations` script. In addition, the {RR_jackknife_counts_XY} quantities are the :math:`RR_{aA}^{XY}` pair counts which can be specified to avoid recomputation. These have been previously output by the :doc:`jackknife-weights` code as ``jackknife_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` (using the correct covariance-matrix binning) hence can be used here for a significant speed boost. The :math:`RR_{aA}^{XY}` pair counts must be normalized by the squared full-survey summed weights :math:`(\sum_i w_i)^2` - this is done automatically in the preceding script.


**Output Files**

This script creates ASCII files for each output correlation function, of the form ``xi_jack_n{N}_m{M}_{INDEX}.dat`` for N radial bins, M angular bins and INDEX specifying the correlation function type (11 = autocorrelation of field 1 (default), 12 = cross-correlation of fields 1 and 2, 22 = autocorrelation of field 2). **NB**: These have a different file format to the non-jackknife correlation functions. The first and second lines of the ``.dat`` file list the radial and angular bin centers, but each succeeding line gives the entire correlation function estimate for a given jackknife. The rows indicate the jackknife and the columns specify the collapsed bin, using the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins. 

These files are read automatically by the :ref:`post-processing-multi` code.
