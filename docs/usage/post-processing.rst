Post-Processing & Reconstruction
=================================

These scripts post-process the single- or multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter(s), :math:`alpha_i`, from data derived covariance matrices (from individual jackknife correlation functions computed in the :ref:`jackknife-correlations` script). A variety of analysis products are output as an ``.npz`` file, as described below.

.. _post-processing-single:

Single-Field Reconstruction
------------------------------

This reconstructs output covariance matrices for a single field. Before running this script, covariance matrix estimates must be produced using the :doc:`main-code` code. The shot-noise rescaling parameters are computed by comparing the theoretical jackknife covariance matrix :math:`\hat{C}^{J}_{ab}(\alpha)` with that computed from the data itself, using individual unrestricted jackknife estimates :math:`\hat{\xi}^J_{aA}`. We define the data jackknife covariance matrix as :math:`C^{\mathrm{data}}_{ab} = \sum_A w_{aA}w_{bA}\left(\hat\xi^J_{aA} - \bar{\xi}_a\right)\left(\hat\xi^J_{bA}-\bar\xi_b\right) / \left(1-\sum_B w_{aB} w_{bB}\right)`, where :math:`\bar\xi_a` is the mean correlation function in bin :math:`a`. We compute :math:`\alpha` via minimizing the likelihood function :math:`-\log\mathcal{L}_1(\alpha) = \mathrm{trace}(\Psi^J(\alpha)C^\mathrm{data}) - \log\mathrm{det}\Psi^J(\alpha)+\mathrm{const}.` using the (bias-corrected) precision matrix :math:`\Psi(\alpha)`.

**Usage**::
    
    python python/post-process.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}

**Input Parameters**

- {XI_JACKKNIFE_FILE}: Input ASCII file containing the correlation function estimates :math:`\xi^J_A(r,\mu)` for each jackknife region, as created by the :ref:`jackknife-correlations` script. This has the format specified in :ref:`file-inputs`.
- {WEIGHTS_DIR}: Directory containing the jackknife weights and pair counts, as created by the :doc:`jackknife-weights` script. This must contain ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` using the relevant covariance matrix binning scheme.
- {COVARIANCE_DIR}: Directory containing the covariance matrix ``.txt`` files output by the :doc:`main-code` C++ code. This directory should contain the subdirectories ``CovMatricesAll`` and ``CovMatricesJack`` containing the relevant analysis products.
- {N_MU_BINS}: Number of angular (:math:`\mu`) bins used in the analysis.
- {N_SUBSAMPLES}: Number of individual matrix subsamples computed in the C++ code. This is the :math:`N_\mathrm{loops}` parameter used in the :ref:`main-code` code. Individual matrix estimates are used to remove quadratic bias from the precision matrix estimates and compute the effective number of degrees of freedom :math:`N_\mathrm{eff}`.
- {OUTPUT_DIR}: Directroy in which to save the analysis products. This will be created if not present.

**Output**

This script creates a single compressed Python file ``Rescaled_Covariance_Matrices_n{N}_m{M}_j{J}.npz`` as an output in the given output directory. All matrices follow the collapsed bin indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins and have dimension :math:`n_\mathrm{bins}\times n_\mathrm{bins}` for a total of :math:`n_\mathrm{bins}` bins. Precision matrices are computed using the quadratic bias elimination method of `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_. All matrices are output using the optimal shot-noise rescaling parameter.  

The output file has the following entries:

- :attr:`shot_noise_rescaling` (Float): Optimal value of the shot-noise rescaling parameter, :math:`\alpha^*`, from the :math:`\mathcal{L}_1` maximization. 
- :attr:`jackknife_theory_covariance` (np.ndarray): Theoretical jackknife covariance matrix estimate :math:`\hat{C}^J_{ab}(\alpha^*)`.
- :attr:`full_theory_covariance` (np.ndarray): Theoretical full covariance matrix estimate :math:`\hat{C}_{ab}(\alpha^*)`.
- :attr:`jackknife_data_covariance` (np.ndarray): Data-derived jackknife covariance matrix :math:`\hat{C}^\mathrm{data}_{ab}`, computed from the individual unrestricted jackknife correlation function estimates.
- :attr:`jackknife_theory_precision` (np.ndarray): Associated precision matrix to the theoretical jackknife covariance matrix estimate, :math:`\Psi_{ab}^J(\alpha^*)`. 
- :attr:`full_theory_precision` (np.ndarray): Associated precision matrix to the theoretical full covariance matrix estimate, :math:`\Psi_{ab}(\alpha^*)`.
- :attr:`individual_theory_covariances` (list): List of individual (and independent) full theoretical covariance matrix estimates. These are used to compute :math:`\tilde{D}_{ab}` and comprise N_SUBSAMPLES estimates.
- :attr:`full_theory_D_matrix` (np.ndarray): Quadratic bias correction :math:`\tilde{D}_{ab}` matrix for the full theoretical covariance matrix, as described in `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_.
- :attr:`N_eff` (Float): Effective number of mocks in the output full covariance matrix, :math:`N_\mathrm{eff}`, computed from :math:`\tilde{D}_{ab}`.


.. _post-processinng-multi:

Multi-Field Reconstruction
-----------------------------

.. todo:: include the multi-field reconstruction codes + estimate of the jackknife data covariance.




Here, we compute the weights assigned to each jackknife region for each bin. This is done using the `Corrfunc <https://corrfunc.readthedocs.io>`_ of Sinha & Garrison to compute the weights :math:`w_{aA}^{XY} = RR_{aA}^{XY} / \sum_B RR_{aB}^{XY}` for bin :math:`a`, jackknife :math:`A` and fields :math:`X` and :math:`Y`. 

Two codes are supplied; one using a single set of tracer particles and the other with two input sets, for computation of cross-covariance matrices. These are in the ``python/`` directory. This must be run before the main C++ code.

Usage
~~~~~~~
For a single field analysis::

    python python/jackknife_weights.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} OUTPUT_DIR}

For an analysis using two distinct fields::

    python python/jackknife_weights_cross.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR}
    
**NB**: The two field script computes all three combinations of weights between the two random fields, thus has a runtime :math:`\sim` 3 times that of ``jackknife_weights.py``. Running these together in one script ensures that we have the same number of jackknives for all fields. Also, the two fields must be distinct, else there are issues with double counting. 

.. todo:: check RascalC read-in procedure with all weights 
    
**Input Parameters**

- {RANDOM_PARTICLE_FILE}, {RANDOM_PARTICLE_FILE_1}, {RANDOM_PARTICLE_FILE_2}: Input ASCII file containing random particle positions and jackknife numbers in {x,y,z,weight,jackknife_ID} format, such as that created with the :doc:`pre-processing` scripts. This should be in ``.csv``, ``.txt`` or ``.dat`` format with space-separated columns.
- {BIN_FILE}: ASCII file specifying the radial bins, as described in :ref:`file-inputs`. This can be user-defined or created by the :ref:`write-binning-file` scripts.
- {MU_MAX}: Maximum :math:`\mu = \cos\theta` used in the angular binning.
- {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu]`.
- {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
- {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
- {OUTPUT_DIR}: Directory in which to house the jackknife weights and pair counts. This will be created if not in existence.


**Notes**:

- This is a very CPU intensive computation since we must compute pair counts between every pair of random particles. The process can be expedited using multiple CPU cores or a reduced number of random particles (e.g. via the :ref:`particle-subset` script).
- For two sets of input particles, three sets of weights must be computed for the three possible pairs of two distinct fields, hence the computation time increases by a factor of three.

**Note on Periodicity**

The code can be run for datasets created with either periodic or non-periodic boundary conditions. Periodic boundary conditions are often found in cosmological simlulations. If periodic, the pair-separation angle :math:`\theta` (used in :math:`\mu=\cos\theta`) is measured from the :math:`z` axis, else it is measured from the radial direction. If periodic data is used, the C++ code **must** be compiled with the -DPERIODIC flag.

Output files
~~~~~~~~~~~~~

This code creates ASCII files containing the jackknife weights for each bin, the RR pair counts for each bin in each jackknife and the summed RR pair counts in each bin. The output files have the format ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat``, ``jackknife_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` respectively N and M specify the number of radial and angular bins respectively and J gives the number of non-empty jackknives. INDEX specifies which fields are being used i.e. INDEX = 12 implies the :math:`w_{aA}^{12}`, :math:`RR_{aA}^{12}` and :math:`RR_a^{12}` quantities.

The binned pair counts is a list of weighted pair counts for each bin, summed over all jackknife regions, in the form :math:`RR_a^{J,XY} = \sum_B RR_{aB}^{XY}`, with each bin on a separate row. The jackknife pair counts and jackknife weights files list the quantities :math:`RR_{aA}^{XY}` and :math:`w_{aA}^{XY}` for each bin and jackknife region respectively. The :math:`j`-th row contains the (tab-separated) quantities for each bin using the :math:`j`-th jackknife. The first value in each row is the jackknife number, and the bins are ordered using the collapsed binning :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins.  
