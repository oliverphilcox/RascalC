Fundamentals
============

Introduction
------------

Here we strive to provide a practical summary.
We also recommend reading a theoretical overview of the RascalC methodology in Sections 2.1 and 3 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_, or Section 2.2 of `Rashkovetskyi 2025 dissertation <https://rashkovetsky.im/files/Dissertation_Rashkovetskyi_2025.pdf>`_.
Many of the more technical (algorithm implementation) details omitted there are described in Sections 3 and 4 of `Philcox et al 2020 <https://arxiv.org/abs/1904.11070>`_.

RascalC constructs a one-parametric covariance matrix model (in the single-tracer case) from an arbitrary (e.g., empirical or theoretical) 2-point correlation function (2PCF):

.. math:: C(\alpha_{\rm SN}) = C_4 + C_3 \alpha_{\rm SN} + C_2 \alpha_{\rm SN}^2
    :label: cov_model

where :math:`C_4, C_3` and :math:`C_2` are 4-, 3- and 2-point terms respectively.
The model parameter :math:`\alpha_{\rm SN}` scales density, or equivalently the amount of shot noise, and is accordingly called **shot-noise rescaling**.
We recommend calibrating this shot-noise rescaling parameter using a reference covariance from jackknife resampling of the data (:ref:`pipeline_jack`), or variations in a sample of mocks (:ref:`pipeline_mock`).

The model is similar to the Gauss-Poisson approximation in Fourier space (see e.g. `Grieb et al 2016 <https://arxiv.org/abs/1509.04293>`_) allowing the shot noise to be a free parameter.
Working in configuration space allows to naturally account for survey geometry and selection, including the variation of the expected density :math:`\bar n` across the survey volume, by sampling points from the random catalog.

**This Python library currently only implements the covariance estimators for 2PCF**.
Covariances for 3PCF (from later parts of `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_) and configuration-space power spectrum estimators (`Philcox & Eisenstein 2020 <https://arxiv.org/abs/1912.01010>`_) are not supported (yet).
To use them, one would probably need to dig in the historic C++ code.

.. _pipeline_basic:

Basic/minimal pipeline
^^^^^^^^^^^^^^^^^^^^^^

Use with caution!
This pipeline does not include a reference to determine the shot-noise rescaling parameter.
The default value for the shot-noise rescaling is 1, but we can not recommend using it in real applications.
Refer to :ref:`pipeline_jack` or :ref:`pipeline_mock` **after looking at the** :ref:`pipeline_basic_fig` **and reading the** :ref:`general_usage` **after it**.
However, we can suggest at least two situations in which this variant is advisable:

    - Two-tracer covariance when the shot-noise rescaling values for each tracer are known (or will be known soon) from single-tracer jackknife computations. Currently, building a two-tracer jackknife model seems excessive.
    - You can also use the basic pipeline to re-run the post-processing with a mock 2PCF sample later. However, substituting jackknife can not be deferred like this because the jackknife model crucially depends on the jackknife RR counts.

.. digraph:: pipeline_basic
    :name: pipeline_basic_fig
    :caption: Basic/minimal pipeline flowchart

    "Random catalog", "Data catalog" [style=filled, fillcolor=red];
    "RR counts", "Full 2PCF", "Shot-noise rescaling" [shape=egg, style=filled, fillcolor=orange];
    "RascalC", "Full covariance model" [shape=box, style=filled, fillcolor=yellow];
    "RascalC" [fontname="Courier New"];
    "Random catalog" -> {"RascalC" "RR counts" "Full 2PCF"};
    "Data catalog" -> "Full 2PCF";
    {"RR counts" "Full 2PCF"} -> "RascalC" -> "Full covariance model";
    {"Full covariance model" "Shot-noise rescaling"} -> "Full (final) covariance";
    "Full (final) covariance" [style=filled, fillcolor=green];

.. _general_usage:

General guidance
^^^^^^^^^^^^^^^^

General practical usage remarks for the Python wrapper function, :func:`RascalC.run_cov`:

- RR counts and 2PCF can and should all be estimated at the same time using the ``pycorr`` `library for 2-point correlation function estimation <https://github.com/cosmodesi/pycorr>`_; this is an external step to :func:`RascalC.run_cov`.

    - Use ``s_mu`` mode in ``pycorr``, other counting modes are not supported by ``RascalC``.
    - Use the Landy-Szalay 2PCF estimator, or the natural 2PCF estimator with analytical RR counts in periodic boxes. Alternative estimators have a different (typically higher) variance and are not supported by ``RascalC``.
    - The necessary pair counts can be computed on GPU, whereas ``RascalC`` can only use CPU (currently).
    - Even if you need to use CPU, **you should run the counts in a separate, independent process from the one calling** :func:`RascalC.run_cov`, **because both should be parallelized and they are known to interfere with each other's efficiency**.
- Choose a binning ``mode`` for the covariance:

    - ``s_mu`` mode for angular bins (uniform in :math:`0 \le \left| \mu \right| \le 1`, where :math:`\mu \equiv \cos \theta`; the mode was implemented in `Philcox et al. 2020 <https://arxiv.org/abs/1904.11070>`_ and used in `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_);
    - ``legendre_projected`` mode for Legendre multipole moments in separation (radial) bins, corresponding to the ``pycorr`` multipole estimation via projection from angular bins (introduced and validated in `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_).
    - ``legendre_accumulated`` mode for Legendre multipole moments accumulated at pair-counting, using a survey correction function with realistic survey geometry (introduced in `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_). It is simpler with periodic cubic boxes, but this mode is not compatible with jackknives.
- Load the ``pycorr`` 2PCF estimator computed beforehand, **cut and/or rebin it to radial (separation) bins desired for the covariance** (e.g., 4 Mpc/h wide from 20 to 200 Mpc/h) and pass it through the ``pycorr_allcounts_11`` argument.

    - In ``s_mu`` mode, you should also rebin angularly to the desired number of angular bins, barring wrapping of :math:`-1 \le \mu < 0`, as explained next:
    
        - It is recommended that you leave the counts in :math:`-1 \le \mu \le 1` bins (for potential error correction), the code will wrap them to :math:`0 \le \left| \mu \right| \le 1` automatically, halving the number of angular bins.
    - In Legendre modes, you can leave the angular (:math:`\mu`) bins as they are.
- You also need to provide input clustering in a form of 2PCF table via the ``xi_table_11`` argument. You can use the pre-computed and loaded ``pycorr`` 2PCF estimator again, but you might want to rebin it differently from the previous case. It is usually advisable to have ``xi_table_11`` in finer radial bins than ``pycorr_allcounts_11``, but the angular (:math:`\mu`) should not be too fine to avoid noisiness.
- Random positions are another necessary input as the ``randoms_positions1`` argument.

    - The number of randoms for ``RascalC`` does not have to be the same as for pair counting and 2PCF estimation (except when you disable ``normalize_wcounts``). It should be high enough to provide a good representation of survey geometry, but not too high to keep the run time reasonable.
    - For 2PCF covariance after standard BAO reconstruction, provide the **shifted** randoms (`Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_, `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_; the input 2PCF conversion in that case will be applied automatically to a ``pycorr`` estimator).
    - For a periodic cubic box (without reconstruction), you will need to generate uniform random positions yourself.
- You must also pass the weights for the randoms through the ``randoms_weights1`` argument. These must match what you used for pair counting and 2PCF estimation with ``pycorr``.

    - If you did not use weights in ``pycorr``, you should pass an array containing 1 for each random point as ``randoms_weights1``.
- Set the output directory, ``out_dir``. We highly recommended a different output directory for each run. This directory will contain all information necessary for post-processing, and a complimentary log file (``log.txt``).
- Set the temporary directory, ``tmp_dir``. Bear in mind that it will need to temporarily contain the random catalog(s) (positions, weights, and jackknife regions if applicable) in text format. This directory can be deleted after the run, and the code normally strives to leave it in its original state. We highly recommended a different temporary directory for each run.
- Set the number of threads via ``nthread``.

    - We recommend trying different options before massive computations.
    - At NERSC Perlmutter, the best option seems to be ``nthread=64`` on half the CPU node (``shared`` queue, requesting 1 node and 128 SLURM cores, which are hyperthreads, and correspond to 64 physical cores).
- Set ``N2``, ``N3`` and ``N4``. These are importance sampling settings: ``RascalC`` will try to take ``N2`` secondary points per each random point provided, ``N3`` tertiary points per each secondary point and ``N4`` quaternary points per each tertiary point. Common values have been ``N2=5, N3=10, N4=20``, but you may need to adjust them (see :ref:`improving_convergence`).
- Set the number of integration loops ``n_loops`` for covariance matrix terms evaluation. It should be divisible by ``nthread``. 1024 might be a nice starting value, but you may need to adjust it (see :ref:`improving_convergence`).

    - The runtime roughly (not exactly) scales as the number of quads per the number of threads, ``N_randoms * N2 * N3 * N4 * n_loops / nthread`` in single-tracer mode. For reference, on NERSC Perlmutter CPU the code processed about 27 millions (``2.7e7``) quads per second per thread (in December 2024).
- Set the number of loops per sample, ``loops_per_sample``. This sets the amount of auxiliary output used almost exclusively for :ref:`quality_control`. ``loops_per_sample`` needs to be a divider of ``n_loops``, and we recommend keeping ``n_loops / loops_per_sample`` (the number of output subsamples) roughly between 10 and 30. Smaller values may require you to wait too long before there is any usable output, or give insufficient information for :ref:`quality_control`. Larger values can lead to too much output.
- To compute a full two-tracer covariance, you need to also provide all of the following:

    - cross-counts as ``pycorr_allcounts_12`` and second tracer auto-counts as ``pycorr_allcounts_22`` (rebinned in the same way as ``pycorr_allcounts_11``);
    - cross-correlation function as ``xi_table_12`` and the second tracer auto-correlation function as ``xi_table_22`` (in the same format as ``xi_table_11``);
    - second tracer random points positions (``randoms_positions2``) and weights (``randoms_weights2``).
- ``RascalC`` in the flowcharts refers to the most computationally intensive steps (implemented in C++), at which the coefficients for the covariance matrix models are evaluated. These coefficients are saved in a ``Raw_Covariance_Matrices*.npz`` file in the chosen output directory.
- Basic/minimal **post-processing** involves substituting a fixed shot-noise rescaling value (or two values in case of two tracers) into the full covariance model to obtain the final covariance. These operations normally are invoked at the end of :func:`RascalC.run_cov`, but they can also be performed separately using :func:`RascalC.post_process_auto`. The results are saved in a ``Rescaled_Covariance_Matrices*.npz`` file in the chosen output directory.

    - External shot-noise rescaling value(s) (e.g. from jackknife or mock calibration) can be supplied via the ``shot_noise_rescaling1`` parameter (for the first tracer; and ``shot_noise_rescaling2`` for the second tracer if applicable).
    
        - Usage of covariance matrices obtained with the default shot-noise rescaling value(s) of 1 is generally not recommended.
        - These ``shot_noise_rescaling`` parameters set the fixed shot-noise rescaling values only for the :ref:`pipeline_basic`. They are ignored in :ref:`pipeline_jack` and :ref:`pipeline_mock`, where the shot-noise rescaling value is instead determined based on a reference covariance.

After the run (execution of :func:`RascalC.run_cov`):

- If it finished normally with no errors, congratulations!
- If it terminated early and/or the job timed out, it is often worth invoking the post-processing with :func:`RascalC.post_process_auto` to see whether sufficiently good results were saved.
- If there was a convergence-related error or warning, please refer to :ref:`improving_convergence`.
- If anything is unclear, please contact the developer.

In any case, take a look at :ref:`quality_control` after the run.
To work with the final results more conveniently, we recommend seeing :ref:`load_export_final_cov`.

After reading to this point, please refer to :ref:`pipeline_jack` or :ref:`pipeline_mock` to see which fits your needs better.

.. _pipeline_jack:

Jackknife pipeline
^^^^^^^^^^^^^^^^^^

The jackknife pipeline allows to obtain the covariance matrix using only the observational data, without simulations (or using only a single simulated realization).
It has been tested most thoroughly (see e.g. `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_), and is showcased in most :doc:`tutorials-examples`, especially :root:`tutorial.ipynb`.

.. digraph:: pipeline_jack
    :name: pipeline_jack_fig
    :caption: Jackknife pipeline flowchart as in Figure 1 from `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_

    "Random catalog", "Data catalog" [style=filled, fillcolor=red];
    "RR counts (full and jackknife)", "Full 2PCF", "Jackknife 2PCF" [shape=egg, style=filled, fillcolor=orange];
    "RascalC", "Full covariance model", "Jackknife covariance model", "Jackknife covariance", "Best-fit shot-noise rescaling" [shape=box, style=filled, fillcolor=yellow];
    "RascalC" [fontname="Courier New"];
    "Random catalog" -> {"RascalC" "RR counts (full and jackknife)" "Full 2PCF" "Jackknife 2PCF"};
    "Data catalog" -> {"Full 2PCF" "Jackknife 2PCF"};
    {"RR counts (full and jackknife)" "Full 2PCF"} -> "RascalC" -> {"Full covariance model" "Jackknife covariance model"};
    "Jackknife 2PCF" -> "Jackknife covariance";
    {"Jackknife covariance model" "Jackknife covariance"} -> "Best-fit shot-noise rescaling";
    {"Full covariance model" "Best-fit shot-noise rescaling"} -> "Full (final) covariance";
    "Full (final) covariance" [style=filled, fillcolor=green];

Practical remarks particular to the jackknife pipeline with :func:`RascalC.run_cov` in addition to the :ref:`general_usage`:

- Jackknife and full RR counts and 2PCF can and should all be estimated at the same time using the ``pycorr`` `library for 2-point correlation function estimation <https://github.com/cosmodesi/pycorr>`_.

    - Remember that **you should run the counts in a separate, independent process from the one calling** :func:`RascalC.run_cov`, **because both should be parallelized and they are known to interfere with each other's efficiency**.
- The jackknife 2PCF will be loaded from the ``pycorr_allcounts_11`` argument (rebinned as explained above).
- Assign the jackknife regions to the random points (``randoms_positions1``) in the same way as for 2PCF and pair counts, and pass the assignment results (jackknife region number for each random point) through the ``randoms_samples1`` argument.

    - Technically, passing the non-``None`` ``randoms_samples1`` argument switches on the jackknife functionality in :func:`RascalC.run_cov`.
- The code produces a separate model for the jackknife covariance, which takes into account correlations between the jackknife regions.
- Jackknife **post-processing** involves fitting the jackknife covariance model to the data jackknife covariance to find the optimal shot-noise rescaling and substituting that value into the full covariance model to obtain the final covariance. These operations normally are invoked at the end of :func:`RascalC.run_cov`, but they can also be performed separately using :func:`RascalC.post_process_auto`. The results are saved in a ``Rescaled_Covariance_Matrices*Jackknife*.npz`` file in the chosen output directory.

Take a look at :ref:`quality_control` after the run.
To work with the final results more conveniently, we recommend seeing :ref:`load_export_final_cov`.

.. _pipeline_mock:

Mock pipeline
^^^^^^^^^^^^^

Tuning the shot-noise rescaling on mocks was the original method in `O'Connell et al. 2016 <https://arxiv.org/abs/1510.01740>`_; it does not require such a large number of realizations as the direct sample covariance estimation from mocks.
This can be seen as a theory-based template smoothing of the mock sample covariance (where the templates are the full ``RascalC`` covariance model terms), reducing the noise.
However, since `O'Connell & Eisenstein 2018 <https://arxiv.org/abs/1808.05978>`_ it has been largely superseded by the idea of using jackknives, which eliminated the need for mocks except for an occasional validation.
Accordingly, there are few good usage examples, but we are working on this.

.. digraph:: pipeline_mock
    :name: pipeline_mock_fig
    :caption: Mock pipeline flowchart as in Figure 2 from `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_

    "Random catalog", "Data catalog", "Mock catalogs" [style=filled, fillcolor=red];
    "RR counts", "Full 2PCF", "Mock 2PCFs" [shape=egg, style=filled, fillcolor=orange];
    "RascalC", "Full covariance model", "Sample covariance", "Best-fit shot-noise rescaling" [shape=box, style=filled, fillcolor=yellow];
    "RascalC" [fontname="Courier New"];
    "Random catalog" -> {"RascalC" "RR counts" "Full 2PCF"};
    "Data catalog" -> "Full 2PCF";
    {"RR counts" "Full 2PCF"} -> "RascalC" -> "Full covariance model";
    "Mock catalogs" -> "Mock 2PCFs" -> "Sample covariance";
    {"Full covariance model" "Sample covariance"} -> "Best-fit shot-noise rescaling";
    {"Full covariance model" "Best-fit shot-noise rescaling"} -> "Full (final) covariance";
    "Full (final) covariance" [style=filled, fillcolor=green];

Currently, the mock pipeline can only be used by

1. running the :ref:`pipeline_basic`;
2. computing and writing the mock 2PCF sample covariance with e.g.

    - :func:`RascalC.pycorr_utils.sample_cov.sample_cov_from_pycorr_to_file` in ``s_mu`` mode;
    - :func:`RascalC.pycorr_utils.sample_cov_multipoles.sample_cov_multipoles_from_pycorr_to_file` in Legendre multipole modes;
3. invoking the manual and rather tedious post-processing by choosing an appropriate function:

    - :func:`RascalC.post_process.post_process_default_mocks` in ``s_mu`` mode (for a single tracer);

        - :func:`RascalC.post_process.post_process_default_mocks_multi` in ``s_mu`` mode for two tracers;
    - :func:`RascalC.post_process.post_process_legendre_mocks` in Legendre modes (for a single tracer).

The post-processing results will be saved in a ``Rescaled_Covariance_Matrices*Mocks*.npz`` file in the chosen output directory.
Take a look at :ref:`quality_control` after the run.
To work with the final results more conveniently, we recommend seeing :ref:`load_export_final_cov`.

We are working on allowing to pass the mock correlation functions and/or the mock sample covariance to :func:`RascalC.run_cov` and :func:`RascalC.post_process_auto` to make this pipeline easier to use.

.. _quality_control:

General quality control
^^^^^^^^^^^^^^^^^^^^^^^

The convergence checks mostly follow Section 6.1 of `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_.

1. The strictest criterion is that the final covariance matrix should be positive definite. If this condition is violated, the Python code raises an exception, which should be easy to notice.
2. Next, there is the eigenvalue test, which produces warnings if failed:

    - In the original (stronger) version (Equation (4.5) in `Philcox et al 2020 <https://arxiv.org/abs/1904.11070>`_), the minimal eigenvalue of the 4-point covariance term :math:`C_4` should be larger than minus the minimal eigenvalue of the 2-point term :math:`C_2`, both in :eq:`cov_model`.

        - Shot-noise rescaling values smaller than 1 (which are quite common, e.g. in `Rashkovetskyi et al 2025 <https://arxiv.org/abs/2404.03007>`_) make this criterion stricter, because they scale the 2-point term down by :math:`\alpha_{\rm SN}^2`. So the code now repeats the test is the optimal shot-noise rescaling value becomes less than 1.
    - On the other hand, the compared eigenvalues of :math:`C_4` and :math:`C_2` can correspond to quite different separation scales, making the original criterion unnecessarily strict in some cases. This consideration led us to introduce the weaker version, where we compare the eigenvalues of :math:`C_2^{-1/2} C_4 C_2^{-1/2}` with :math:`-1` or :math:`-\alpha_{\rm SN}^2`. Here :math:`C_2^{-1/2}` is the inverse of `the matrix square root <https://en.wikipedia.org/wiki/Square_root_of_a_matrix>`_ (`scipy.linalg.sqrtm <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.sqrtm.html>`_) of the 2-point term, which scales the different parts of the 4-point term matrix more appropriately. (The 2-point term is either diagonal or block-diagonal with small blocks, so taking its matrix square root should be numerically stable.)
3. Finally, there is the extra convergence check (:mod:`RascalC.convergence_check_extra`) performed at the end of :func:`RascalC.run_cov` or :func:`RascalC.post_process_auto` by default.

    - After Section 3.2 of `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_, we recommend focusing on ``R_inv`` (:math:`R_{\rm inv}`) values. There is no universal threshold, but some decent reference values are:

        - <0.6% (``6e-3``) for Early DESI data BGS/LRG with 45 bins (45 radial times 1 angular);
        - <5% (``5e-2``) for DESI DR1/DR2 LRG/ELG, <12% (``1.2e-1``) for ``BGS_BRIGHT-21.5`` (:math:`0.1<z<0.4`) and <20% (``2e-1``) for ``BGS_BRIGHT-21.35`` with 135 bins (45 radial times 3 multipoles). QSO have almost always converged much better, like 0.1-0.2% (``2e-3``) due to higher shot-noise, making the easy 2-point term the dominant one.
        - Values exceeding 1 are quite certainly high.
    - These figures of intrinsic scatter in covariance sums/integrals estimated with importance sampling tend to increase

        - as the number of bins increases (the trend is the same for mocks — see e.g. Equation (3.12) in `Rashkovetskyi et al 2023 <https://arxiv.org/abs/2306.06320>`_)
        - as the sample density increases and shot-noise decreases (parallel with mocks is less clear, but dense samples are also harder to simulate).

.. _improving_convergence:

Addressing convergence issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these instructions when

- you get the error message ``The full covariance is not positive definite - insufficient convergence``;
- the ``R_inv`` values from the extra convergence check are worryingly high (see above for reference, or reach out if in doubt);
- you get the ``4-point covariance matrix has not converged properly via the weaker eigenvalue test`` warnings, although they seldom appear without one of the previous two issues;

    - you can probably ignore the warning(s) about ``the stronger eigenvalue test`` if they appear alone.

First, you can re-run post-processing with alternative options using :func:`RascalC.post_process_auto`:

- skip a few bins with smallest separations by passing a single positive integer (their number) via ``skip_s_bins``;

    - also skip a few bins with highest separations by passing a tuple of two positive numbers via ``skip_s_bins``; among them, the first number sets how many bins to skip at the low end, the second — at the high end;
- in Legendre mode, you can also try skipping highest multipoles by passing their number (counting only even multipoles) via ``skip_l``.
- you can also try to discard "unlucky" samples using the ``n_samples`` argument, but this seldom helps and can become confusing.

The above are the fastest options because they only require re-running the post-processing script while re-using the products of the main computation.

Then, consider running the main computation again with :func:`RascalC.run_cov` to generate a smaller covariance matrix by

- using coarser radial (separation) bins;
- in ``s_mu`` mode, using coarser angular (:math:`\mu`) bins;
- in Legendre mode, using fewer multipoles (this alone can be done faster and easier in post-processing, as was just mentioned above, but not in combination with different radial/separation bins).

If changing the covariance binning is undesirable or does not help, you should run longer.
Use a different output directory to prevent confusion and overrides; renaming the old output directory also works.
This can be done in different ways:

- Increase the number of loops (``n_loops``).

    - Occasionally bad convergence is just bad luck, so running again with the same settings, including ``n_loops`` might not be needed. In that case, just do not use the same fixed ``seed``, as that should reproduce the results exactly.
    - If you keep other settings fixed (except ``n_loops``, ``nthread`` and naturally ``seed`` — you can change those more freely), you can also concatenate (combine) samples from different runs into a new, different directory using :func:`RascalC.cat_raw_covariance_matrices` to reach even better convegence (check it by post-processing the new directory with :func:`RascalC:post_process_auto`). However, combining samples does not always improve convergence, and keeping track of different sample combination can be hard.
- Increase ``N2`` — probably not recommended, because the effect is similar to increasing ``n_loops``, but sample combination is no longer an optionl.
- Increase ``N4`` and/or ``N3``. It is probably more sensible than the previous options because we expect the higher-point terms to converge slower. ``N4`` will only affect the 4-point term :math:`C_4`; ``N3`` also affects the 3-point term :math:`C_3`.

The convergence issues can be persistent.
Please do not hesitate to reach out.

The main computation wrapper
----------------------------

This Python wrapper of ``RascalC`` is heavily interfaced with ``pycorr`` `library for 2-point correlation function estimation <https://github.com/cosmodesi/pycorr>`_.
Many of the arguments are intentionally similar to ``pycorr.TwoPointCorrelationFunction`` `high-level interface <https://py2pcf.readthedocs.io/en/latest/api/api.html#pycorr.correlation_function.TwoPointCorrelationFunction>`_.

Please bear with the long description; you can pay less attention to settings labeled optional in the beginning.

.. autofunction:: RascalC.run_cov


Post-processing
---------------

A suitable post-processing routine is invoked at the end of the main wrapper function (:func:`RascalC.run_cov`), so in many circumstances you may not need to run it separately.
However, this automated but customizable post-processing routine is useful for timed-out runs, switching the mode, testing different cuts and/or output combinations in cases of insufficient convergence, etc.

.. autofunction:: RascalC.post_process_auto

.. _load_export_final_cov:

Loading and exporting the final covariance matrices
---------------------------------------------------

Perform the :ref:`quality_control` first.

.. automodule:: RascalC.cov_utils
    :members:
    :member-order: bysource