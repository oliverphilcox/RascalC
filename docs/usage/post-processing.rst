Post-Processing & Reconstruction
=================================

These scripts post-process the single- or multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameter(s), :math:`\alpha_i`, from data derived covariance matrices (from individual jackknife correlation functions computed in the :ref:`jackknife-correlations` script). A variety of analysis products are output as an ``.npz`` file, as described below.


.. todo:: test eigenvalue test for convergence in both scripts

.. todo:: test full multi-field post-processing scripts

.. _post-processing-single:


Single-Field Reconstruction
------------------------------

This reconstructs output covariance matrices for a single field. Before running this script, covariance matrix estimates must be produced using the :doc:`main-code` code. The shot-noise rescaling parameters are computed by comparing the theoretical jackknife covariance matrix :math:`\hat{C}^{J}_{ab}(\alpha)` with that computed from the data itself, using individual unrestricted jackknife estimates :math:`\hat{\xi}^J_{aA}`. We define the data jackknife covariance matrix as :math:`C^{\mathrm{data}}_{ab} = \sum_A w_{aA}w_{bA}\left(\hat\xi^J_{aA} - \bar{\xi}_a\right)\left(\hat\xi^J_{bA}-\bar\xi_b\right) / \left(1-\sum_B w_{aB} w_{bB}\right)`, where :math:`\bar\xi_a` is the mean correlation function in bin :math:`a`. We compute :math:`\alpha` via minimizing the likelihood function :math:`-\log\mathcal{L}_1(\alpha) = \mathrm{trace}(\Psi^J(\alpha)C^\mathrm{data}) - \log\mathrm{det}\Psi^J(\alpha)+\mathrm{const}.` using the (bias-corrected) precision matrix :math:`\Psi(\alpha)`.

**NB**: Before full computation of precision matrices and shot-noise rescaling, we check that the matrices have converged to a sufficient extent to allow convergence. As described in the RascalC paper, we require :math:`\text{min\,eig}(C_4) \geq - \text{min\,eig}(C_2)` to allow inversion. If this condition is met (due to too small sampling time, i.e. too low :math:`N_i` and/or :math:`N_\mathrm{loops}` values), the code exits.

**Usage**::
    
    python python/post-process.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}

**Input Parameters**

- {XI_JACKKNIFE_FILE}: Input ASCII file containing the correlation function estimates :math:`\xi^J_A(r,\mu)` for each jackknife region, as created by the :ref:`jackknife-correlations` script. This has the format specified in :ref:`file-inputs`.
- {WEIGHTS_DIR}: Directory containing the jackknife weights and pair counts, as created by the :doc:`jackknife-weights` script. This must contain ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` using the relevant covariance matrix binning scheme.
- {COVARIANCE_DIR}: Directory containing the covariance matrix ``.txt`` files output by the :doc:`main-code` C++ code. This directory should contain the subdirectories ``CovMatricesAll`` and ``CovMatricesJack`` containing the relevant analysis products.
- {N_MU_BINS}: Number of angular (:math:`\mu`) bins used in the analysis.
- {N_SUBSAMPLES}: Number of individual matrix subsamples computed in the C++ code. This is the :math:`N_\mathrm{loops}` parameter used in the :ref:`main-code` code. Individual matrix estimates are used to remove quadratic bias from the precision matrix estimates and compute the effective number of degrees of freedom :math:`N_\mathrm{eff}`.
- {OUTPUT_DIR}: Directory in which to save the analysis products. This will be created if not present.

**NB**: This script may take several minutes (to hours) to run if :math:`N_\mathrm{loops}` is larger than a few 10s.

**Output**

This script creates a single compressed Python file ``Rescaled_Covariance_Matrices_n{N}_m{M}_j{J}.npz`` as an output in the given output directory. All matrices follow the collapsed bin indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins and have dimension :math:`n_\mathrm{bins}\times n_\mathrm{bins}` for a total of :math:`n_\mathrm{bins}` bins. Precision matrices are computed using the quadratic bias elimination method of `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_. All matrices are output using the optimal shot-noise rescaling parameter. This file may be large (up to 1GB) depending on the values of :math:`n_\mathrm{bins}` and :math:`N_\mathrm{loops}` used.

The output file has the following entries:

- :attr:`shot_noise_rescaling` (Float): Optimal value of the shot-noise rescaling parameter, :math:`\alpha^*`, from the :math:`\mathcal{L}_1` maximization. 
- :attr:`jackknife_theory_covariance` (np.ndarray): Theoretical jackknife covariance matrix estimate :math:`\hat{C}^J_{ab}(\alpha^*)`.
- :attr:`full_theory_covariance` (np.ndarray): Theoretical full covariance matrix estimate :math:`\hat{C}_{ab}(\alpha^*)`.
- :attr:`jackknife_data_covariance` (np.ndarray): Data-derived jackknife covariance matrix :math:`\hat{C}^{J,\mathrm{data}}_{ab}`, computed from the individual unrestricted jackknife correlation function estimates.
- :attr:`jackknife_theory_precision` (np.ndarray): Associated precision matrix to the theoretical jackknife covariance matrix estimate, :math:`\Psi_{ab}^J(\alpha^*)`. 
- :attr:`full_theory_precision` (np.ndarray): Associated precision matrix to the theoretical full covariance matrix estimate, :math:`\Psi_{ab}(\alpha^*)`.
- :attr:`individual_theory_covariances` (list): List of individual (and independent) full theoretical covariance matrix estimates. These are used to compute :math:`\tilde{D}_{ab}` and comprise N_SUBSAMPLES estimates.
- :attr:`full_theory_D_matrix` (np.ndarray): Quadratic bias correction :math:`\tilde{D}_{ab}` matrix for the full theoretical covariance matrix, as described in `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_.
- :attr:`N_eff` (Float): Effective number of mocks in the output full covariance matrix, :math:`N_\mathrm{eff}`, computed from :math:`\tilde{D}_{ab}`.


.. _post-processing-multi:


Multi-Field Reconstruction
-----------------------------

Analogous to the above, this code performs reconstruction of the covariance matrices, :math:`C_{ab}^{XY,ZW}` for two field cases, using the relevant jackknife correlation functions :math:`\xi^{J,XY}_{aA}` and covariance matrix components. Here, we estimate the shot-noise parameters :math:`\alpha_1` and :math:`\alpha_2` purely from the (11,11) and (22,22) autocovariance matrices, as these give the strongest constraints. In this case, the code will exit if the :math:`C_4^{11,11}` and/or :math:`C_4^{22,22}` are not sufficiently converged, (checking these matrices since :math:`C^{11,11}` and :math:`C^{22,22}` must be inverted to compute :math:`\alpha_1` and :math:`\alpha_2`).

**Usage**::
 
    python python/post_process_multi.py {XI_JACKKNIFE_FILE_11} {XI_JACKKNIFE_FILE_12} {XI_JACKKNIFE_FILE_22} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}

Input parameters are as before, with the addition of :math:`\xi^{J,12}_{aA}` and :math:`\xi^{J,22}_{aA}` files.

**Output**

As above, we create a single compressed Python file for the output analysis products, now labelled ``Rescaled_Multi_Field_Covariance_Matrices_n{N}_m{M}_j{J}.npz``, which contains output matrices for all combinations of the two fields. This could be a large file. This file has the same columns as the single field case, but now :attr:`shot_noise_rescaling` becomes a length-2 array :math:`(\alpha_1^*,\alpha_2^*)`. All other products are are arrays of matrices (shape :math:`2\times2\times2\times2\times n_\mathrm{bins} \times n_\mathrm{bins}`) which are specified by 4 input parameters, corresponding to the desired X, Y, Z, W fields in :math:`C^{XY,ZW}`. This uses Pythonic indexing from 0 to label the input fields. For example, we can access the :math:`\Psi^{11,21}_{ab}` precision matrix by loading the relevant column and specifying the index [0,0,1,0] e.g. to load this matrix we simply use::

    >>> dat=np.load("Rescaled_Multi_Field_Covariance_Matrices_n36_m12_j169.npz") # load the full data file
    >>> full_precision = dat['full_theory_precision'] # load the precision matrix
    >>> psi_1121 = full_precision[0,0,1,0] # specify the (11,21) component
