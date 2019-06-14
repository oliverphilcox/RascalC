Post-Processing & Reconstruction
=================================

These scripts post-process the single- or multi-field integrals computed by the C++ code described in :doc:`main-code`. For the DEFAULT, LEGENDRE and 3PCF modes, the output covariance matrix products (as well as precision matrices and effective number of mocks) are computed and saved to a ``.npz`` file, as described below. This is done for some input shot-noise rescaling parameter :math:`\alpha` (or one for each field, for the multi-field case).

In JACKKNIFE mode, the shot-noise rescaling parameter(s), :math:`\alpha_i`, are computed from data derived covariance matrices (from individual jackknife correlation functions computed in the :ref:`jackknife-correlations` script). As above, a variety of analysis products are output as an ``.npz`` file.

**NB**: Before full computation of precision matrices and shot-noise rescaling, we check that the matrices have converged to a sufficient extent to allow convergence. As described in the RascalC paper, we require :math:`\text{min eig}(C_4) \geq - \text{min eig}(C_2)` to allow inversion. If this condition is met (due to too small sampling time, i.e. too low :math:`\{N_i\}` and/or :math:`N_\mathrm{loops}` values), the code exits. For multi-field cases, we require the (11,11) and (22,22) autocovariance matrices to be sufficiently converged.

**General Input Parameters**

- {COVARIANCE_DIR}: Directory containing the covariance matrix ``.txt`` files output by the :doc:`main-code` C++ code. This directory should contain the subdirectories ``CovMatricesAll`` and ``CovMatricesJack`` containing the relevant analysis products.
- {N_R_BINS}: Number of radial  bins used in the analysis.
- {N_MU_BINS}: Number of angular (:math:`\mu`) bins used in the analysis.
- {MAX_L}: Maximum Legendre multipole used in the analysis.
- {N_SUBSAMPLES}: Number of individual matrix subsamples computed in the C++ code. This is the :math:`N_\mathrm{loops}` parameter used in the :ref:`main-code` code. Individual matrix estimates are used to remove quadratic bias from the precision matrix estimates and compute the effective number of mocks :math:`N_\mathrm{eff}`.
- {OUTPUT_DIR}: Directory in which to save the analysis products. This will be created if not present.
- {SHOT_NOISE_RESCALING}, {SHOT_NOISE_RESCALING_1}, {SHOT_NOISE_RESCALING_2} *(DEFAULT, LEGENDRE and 3PCF modes)*: Shot-noise rescaling parameter to be applied to the output covariance and precision matrices. This also affects the effective number of mocks, :math:`N_\mathrm{eff}`. If unspecified, this is set to unity.


**Output**

The single field scripts create a single compressed Python file ``Rescaled_Covariance_Matrices_{MODE}_n{N}_m{M}[_j{J}].npz`` as an output in the given output directory. All 2PCF (3PCF) matrices follow the collapsed bin indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` (:math:`\mathrm{bin}_\mathrm{collapsed} = \left(\mathrm{bin}_{\mathrm{radial},1}\times n_r + \mathrm{bin}_{\mathrm{radial},2}\right)\times n_\mu + \mathrm{bin}_\mathrm{angular}`) for a total of :math:`n_\mu` angular bins and have dimension :math:`n_\mathrm{bins}\times n_\mathrm{bins}` for a total of :math:`n_\mathrm{bins}` bins. Precision matrices are computed using the quadratic bias elimination method of `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_. All matrices are output using the optimal shot-noise rescaling parameter. This file may be large (up to 1GB) depending on the values of :math:`n_\mathrm{bins}` and :math:`N_\mathrm{loops}` used.

The output file has the following entries:

- :attr:`shot_noise_rescaling` (Float): Optimal value of the shot-noise rescaling parameter, :math:`\alpha^*`, from the :math:`\mathcal{L}_1` maximization. 
- :attr:`full_theory_covariance` (np.ndarray): Theoretical full covariance matrix estimate :math:`\hat{C}_{ab}(\alpha^*)`.
- :attr:`full_theory_precision` (np.ndarray): Associated precision matrix to the theoretical full covariance matrix estimate, :math:`\Psi_{ab}(\alpha^*)`.
- :attr:`individual_theory_covariances` (list): List of individual (and independent) full theoretical covariance matrix estimates. These are used to compute :math:`\tilde{D}_{ab}` and comprise N_SUBSAMPLES estimates.
- :attr:`full_theory_D_matrix` (np.ndarray): Quadratic bias correction :math:`\tilde{D}_{ab}` matrix for the full theoretical covariance matrix, as described in `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_.
- :attr:`N_eff` (Float): Effective number of mocks in the output full covariance matrix, :math:`N_\mathrm{eff}`, computed from :math:`\tilde{D}_{ab}`.

For multi-field cases, we also create a single compressed Python file for the output analysis products, now labelled ``Rescaled_Multi_Field_Covariance_Matrices_{MODE}_n{N}_m{M}[_j{J}].npz``, which contains output matrices for all combinations of the two fields. This could be a large file. This file has the same columns as the single field case, but now :attr:`shot_noise_rescaling` becomes a length-2 array :math:`(\alpha_1^*,\alpha_2^*)`. All other products are are arrays of matrices (shape :math:`2\times2\times2\times2\times n_\mathrm{bins} \times n_\mathrm{bins}`) which are specified by 4 input parameters, corresponding to the desired X, Y, Z, W fields in :math:`C^{XY,ZW}`. This uses Pythonic indexing from 0 to label the input fields. For example, we can access the :math:`\Psi^{11,21}_{ab}` precision matrix by loading the relevant column and specifying the index [0,0,1,0] e.g. to load this matrix we simply use::

    >>> dat=np.load("Rescaled_Multi_Field_Covariance_Matrices_Jackknife_n36_m12_j169.npz") # load the full data file
    >>> full_precision = dat['full_theory_precision'] # load the precision matrix
    >>> psi_1121 = full_precision[0,0,1,0] # specify the (11,21) component


.. _post-processing-general:
    
DEFAULT, LEGENDRE and 3PCF mode reconstruction
-----------------------------------------------

Here we reconstruct the output covariance matrices and associated products, given an input shot-noise rescaling parameter. In 3PCF mode, we do not include the first six-point term, :math:`{}_A^6\mathbf{C}`, as noted in Philcox & Eisenstein (in prep.) since this is expected to be small for a large survey, yet difficult to accurately measure.

**NB**: In 3PCF mode, we require a long integration time for matrix convergence with even a moderate number of bins. If the matrix is not well converged (and invertible) the script will report a precision matrix and effective number of mocks of zero. In many analyses (e.g. Philcox & Eisenstein (in prep.)) the 3PCF covariance is compressed via some metric which improves the convergence. Thus, whilst the full matrix may not be invertible, the compressed version often will be. -

**Usage**

For a single field::

    python python/post_process_default.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]
    python python/post_process_legendre.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]
    python python/post_process_3pcf.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING}]
    
For multiple fields::

    python python/post_process_default_multi.py {COVARIANCE_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2}]
    python python/post_process_legendre_multi.py {COVARIANCE_DIR} {N_R_BINS} {MAX_L} {N_SUBSAMPLES} {OUTPUT_DIR} [{SHOT_NOISE_RESCALING_1} {SHOT_NOISE_RESCALING_2}]


.. _post-processing-jackknife:

JACKKNIFE mode reconstruction
------------------------------

**NB**: This can only be run if the C++ code was run in JACKKNIFE mode for the 2PCF.

This script differs from the above in that we now compute the shot-noise rescaling parameters by comparing the theoretical jackknife covariance matrix :math:`\hat{C}^{J}_{ab}(\alpha)` with that computed from the data itself, using individual unrestricted jackknife estimates :math:`\hat{\xi}^J_{aA}`. We define the data jackknife covariance matrix as :math:`C^{\mathrm{data}}_{ab} = \sum_A w_{aA}w_{bA}\left(\hat\xi^J_{aA} - \bar{\xi}_a\right)\left(\hat\xi^J_{bA}-\bar\xi_b\right) / \left(1-\sum_B w_{aB} w_{bB}\right)`, where :math:`\bar\xi_a` is the mean correlation function in bin :math:`a`. We compute :math:`\alpha` via minimizing the likelihood function :math:`-\log\mathcal{L}_1(\alpha) = \mathrm{trace}(\Psi^J(\alpha)C^\mathrm{data}) - \log\mathrm{det}\Psi^J(\alpha)+\mathrm{const}.` using the (bias-corrected) precision matrix :math:`\Psi^J(\alpha)`. When run for multiple input fields, the (11,11) and (22,22) covariance matrices are used to constrain :math:`\alpha_1` and :math:`\alpha_2` respectively.

**Usage**
    
For a single field::
    
    python python/post-process_jackknife.py {XI_JACKKNIFE_FILE} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}

For multiple fields::

    python python/post_process_jackknife_multi.py {XI_JACKKNIFE_FILE_11} {XI_JACKKNIFE_FILE_12} {XI_JACKKNIFE_FILE_22} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}
    
    
**Additional Jackknife Input Parameters**

- {XI_JACKKNIFE_FILE}, {XI_JACKKNIFE_FILE_11}, {XI_JACKKNIFE_FILE_12}, {XI_JACKKNIFE_FILE_22}: Input ASCII file containing the correlation function estimates :math:`\xi^J_A(r,\mu)` for each jackknife region, as created by the :ref:`jackknife-correlations` script. This has the format specified in :ref:`file-inputs`.
- {WEIGHTS_DIR}: Directory containing the jackknife weights and pair counts, as created by the :doc:`jackknife-weights` script. This must contain ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` using the relevant covariance matrix binning scheme.

**Output**

The output ``.npz`` file contains the following additional columns;

- :attr:`jackknife_theory_covariance` (np.ndarray): Theoretical jackknife covariance matrix estimate :math:`\hat{C}^J_{ab}(\alpha^*)`.
- :attr:`jackknife_data_covariance` (np.ndarray): Data-derived jackknife covariance matrix :math:`\hat{C}^{J,\mathrm{data}}_{ab}`, computed from the individual unrestricted jackknife correlation function estimates.
- :attr:`jackknife_theory_precision` (np.ndarray): Associated precision matrix to the theoretical jackknife covariance matrix estimate, :math:`\Psi_{ab}^J(\alpha^*)`. 


