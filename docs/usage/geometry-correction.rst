Computing Random Counts and Survey Correction Functions
========================================================

Here, we describe the scripts provided to compute RR and RRR random counts, in addition to the *survey correction functions*, defined as the ratio of ideal to true RR (RRR) counts for the 2PCF (3PCF). The RR counts are computed using the `Corrfunc <https://corrfunc.readthedocs.io>`_ code of Sinha & Garrison. We additionally provide functionality to compute multi-field RR and survey-correction functions for the 2PCF. Most of the scripts below are in the ``python/`` directory, and must be run before the main C++ code.

**Note on Periodicity**

The codes below can be run for datasets created with either periodic or non-periodic boundary conditions. Periodic boundary conditions are often found in cosmological simlulations. If periodic, the pair-separation angle :math:`\theta` (used in :math:`\mu=\cos\theta`) is measured from the :math:`z` axis, else it is measured from the radial direction. If periodic data is used, the C++ code **must** be compiled with the -DPERIODIC flag.


.. _RR_counts:

Estimating RR pair counts
--------------------------

This is required to normalize the covariances in DEFAULT mode, and to compute the survey-correction-function in LEGENDRE mode (for **aperiodic** surveys only). In JACKKNIFE mode, the scripts in :doc:`jackknife-weights` should be used instead.

Usage
~~~~~~~
For a single field analysis::

    python python/RR_counts.py {RANDOM_PARTICLE_FILE} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} {NORMED}

For an analysis using two distinct fields::

    python python/RR_counts_multi.py {RANDOM_PARTICLE_FILE_1} {RANDOM_PARTICLE_FILE_2} {BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} {NORMED}

**Input Parameters**

- {RANDOM_PARTICLE_FILE}, {RANDOM_PARTICLE_FILE_1}, {RANDOM_PARTICLE_FILE_2}: Input ASCII file containing random particle positions and jackknife numbers in {x,y,z,weight,[jackknife_ID]} format, such as that created with the :doc:`pre-processing` scripts. This should be in ``.csv``, ``.txt`` or ``.dat`` format with space-separated columns. The 'jackknife_id' column will not be used if included.
- {BIN_FILE}: ASCII file specifying the radial bins, as described in :ref:`file-inputs`. This can be user-defined or created by the :ref:`write-binning-file` scripts.
- {MU_MAX}: Maximum :math:`\mu = \cos\theta` used in the angular binning.
- {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu]`. *NB*: This should be large (preferably :math:`\sim 100`) if the RR counts are being used for LEGENDRE mode survey correction computation. In DEFAULT mode, this should match the output number of angular bins required.
- {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
- {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
- {OUTPUT_DIR}: Directory in which to house the RR counts. This will be created if not in existence.
- {NORMED}: Whether to normalize the output RR counts by the summed squared random weights (0 = no normalization, 1 = normalization). This should be done in LEGENDRE mode (where the RR counts are used to compute the survey correction function) but *not* in the DEFAULT mode.

**Notes**:

- Output files will be saved as ``RR_counts_n{N}_m{M}_{INDEX}.dat`` for N radial and M angular bins. INDEX specifies which fields are used (e.g. INDEX=12 for the field 1 x field 2 pair count). This is saved as a simple list of :math:`N\times M` values using the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}*M + \mathrm{bin}_\mathrm{angular}`.
- This is a very CPU intensive computation since we must compute pair counts between every pair of random particles up to the maximum particle separation. The process can be expedited using multiple CPU cores or a reduced number of random particles (e.g. via the :ref:`particle-subset` script).
- For two sets of input particles, three sets of weights must be computed for the three possible pairs of two distinct fields, hence the computation time increases by a factor of three.


.. _survey_correction_2PCF:

Computing 2PCF Survey Correction Functions
-------------------------------------------

This script computes the survey correction function :math:`\Phi(r_a,\mu)` in each bin :math:`a` and fits it to a smooth model. The output parameters can be fed into the main C++ script. For convenience, we divide the output :math:`\Phi` by :math:`V\overline{(nw)^2}`, for survey volume :math:`V`. For the periodic case, the survey correction function is simply unity everywhere, so it is simple to compute. If aperiodic, we require the survey RR pair counts as an input, either from the :ref:`RR_counts` scripts or elsewhere.

Usage
~~~~~~

For a single field analysis::

    python python/compute_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS}]

For an analysis using two distinct fields::

    python python/compute_correction_function_multi.py {GALAXY_FILE} {GALAXY_FILE_2} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RR_COUNTS_11} {RR_COUNTS_12} {RR_COUNTS_22}]

**Input Parameters**:

See the :ref:`RR_counts` parameters above. Additionally;

- {RR_COUNTS}, {RR_COUNTS_11}, {RR_COUNTS_12} {RR_COUNTS_22} *(Only required if dataset is aperiodic)*: RR counts computed by the :ref:`RR_counts` script, or externally. These should be NOT normalized by the sum of the squared random weights (by setting the NORMED flag above to 0).

**Notes**:

- **NB:** For aperiodic data, this assumes that the weights are FKP weights, such that they can be used to find the random number density at each galaxy position. This is not assumed for periodic data (where number density is constant everywhere).
- Output files are saved as ``BinCorrectionFactor_n{N}_m{M}_{INDEX}.txt`` (aperiodic) or ``BinCorrectionFactor_n{N}_periodic_{INDEX}.txt`` (periodic) for N radial (and M angular) bins, with INDEX specfiying which fields are used. If run for two distinct fields, three correction factors are output, for each non-trivial combination of bins.
- File format is a list of N rows (for N radial bins) with 7 columns specfiying the fitting paramters (as specified in Philcox & Eisenstein 2019). This is automatically read by the main C++ code.
- For a periodic data-set, we output a set of parameters which will lead to the survey correction function being reconstructed as unity everywhere.


.. _RRR_counts:

Estimating RRR triple counts
-----------------------------

We provide a simple code to compute the RRR triple counts, required to find the 3PCF survey correction function for aperiodic survey geometries. (This is not required for periodic surveys). This is a stochastic estimator (similar to the main C++ covariance matrix code) and does *not* exhaustively count all triples of particles. For heavy usage, an exhaustive code such as the :math:`\mathcal{O}(N^2)` algorithm of `Slepian & Eisenstein (2015) <https://arxiv.org/abs/1506.02040>`_ would be more appropriate.

This is a C++ code and is run analogous to the covariance code in :doc:`main-code`. To compile and run use the following:

.. code-block:: bash

    cd triple_counts
    make
    ./triple [OPTIONS]

In the associated Makefile we can specify the ``-DPERIODIC`` and ``-DOPENMP`` flags for periodic and OpenMP parallelized compilation respectively. The command line options (and the majority of the code) are the same as for the main code (see descriptions in :doc:`main-code`), and we must specify at least the input random file, :math:`(r,\mu)` binning strategy and accuracy parameters. *NB*: We must set the minimum :math:`\mu` to be -1 here, as there is no :math:`\mu\leftrightarrow-\mu` symmetry unlike for the 2PCF. In addition, the number of angular bins should be relatively large (:math:`\gtrsim 20`) to reduce bias in later conversion of the correction function into Legendre multipole space.

The code outputs a single text file in the specified output directory, named ``RRR_counts_n{N}_m{M}_full.txt`` which contains the RRR counts, normalized by the summed cubed random particle weights. Each estimate is on a separate line, with the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \left(\mathrm{bin}_{\mathrm{radial},1}\times N + \mathrm{bin}_{\mathrm{radial},2}\right)\times M + \mathrm{bin}_\mathrm{angular}`. These can be read-in by the survey correction function.


.. _survey_correction_3PCF:

Computing 3PCF Survey Correction Functions
-------------------------------------------

This script computes the 3-point survey correction function :math:`\Phi(r_a,r_b,\chi)` in each radial bin pair :math:`a,b`, and outputs the multipoles of :math:`\Phi^{-1}` which are reconstructed by the main C++ code in runtime. For convenience, we divide the output :math:`\Phi` by :math:`6V\overline{(nw)^3}`, for survey volume :math:`V`. For the periodic case, the survey correction function is simply unity everywhere, so it is simple to compute. If aperiodic, we require the (normalized) survey RRR pair counts as an input, either from the :ref:`RRR_counts` scripts or elsewhere.

**NB**: We do not provide functionality for the multi-field 3PCF covariances here.

Usage
~~~~~~
 ::

    python python/compute_3pcf_correction_function.py {GALAXY_FILE} {BIN_FILE} {OUTPUT_DIR} {PERIODIC} [{RRR_COUNTS}]

**Input Parameters**:

See the :ref:`RR_counts` parameters above. Additionally;

- {RRR_COUNTS} *(Only required if dataset is periodic)*: RRR counts computed by the :ref:`RRR_counts` script, or externally. These should be normalized by the sum of the cubed random weights (done natively by the above code).

**Notes**:

- **NB:** For aperiodic data, this assumes that the weights are FKP weights, such that they can be used to find the random number density at each galaxy position. This is not assumed for periodic data (where number density is constant everywhere).
- Output files are saved as ``BinCorrectionFactor3PCF_n{N}_m{M}.txt`` (aperiodic) or ``BinCorrectionFactor3PCF_n{N}_periodic.txt`` (periodic) for N radial (and M angular) bins.
- File format is a list of :math:`N^2` rows (for N radial bins) with 7 columns specfiying the first seven multipoles of the (normalized) inverse :math:`\Phi` multipoles. This is automatically read by the main C++ code.
- For a periodic data-set, we output a set of parameters which will lead to the survey correction function being reconstructed as unity everywhere.
