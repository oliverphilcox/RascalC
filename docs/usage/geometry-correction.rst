Computing Random Counts and the Survey Correction Function
============================================================

.. todo:: Link this into the indexes + getting started etc. Also add 3PCF implementation

Here, we describe the scripts provided to compute RR and RRR random counts, in addition to the *survey correction functions*, defined as the ratio of ideal to true RR (RRR) counts for the 2PCF (3PCF). The RR counts are computed using the `Corrfunc <https://corrfunc.readthedocs.io>`_ code of Sinha & Garrison. We additionally provide functionality to compute multi-field RR and survey-correction functions for the 2PCF. Most of the scripts below are in the ``python/`` directory, and must be run before the main C++ code.

**Note on Periodicity**

The codes below can be run for datasets created with either periodic or non-periodic boundary conditions. Periodic boundary conditions are often found in cosmological simlulations. If periodic, the pair-separation angle :math:`\theta` (used in :math:`\mu=\cos\theta`) is measured from the :math:`z` axis, else it is measured from the radial direction. If periodic data is used, the C++ code **must** be compiled with the -DPERIODIC flag.


.. _RR_counts:

Estimating RR pair counts
--------------------------

This is required to normalize the covariances in DEFAULT mode, and to compute the survey-correction-function in LEGENDRE mode. In JACKKNIFE mode, the scripts in :doc:`jackknife-weights` should be used instead.

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
- {N_MU_BINS}: Number of angular bins used in the range :math:`[0,\mu]`.
- {NTHREADS}: Number of CPU threads to use for pair counting parallelization.
- {PERIODIC}: Whether the input dataset has periodic boundary conditions (0 = non-periodic, 1 = periodic). See note below.
- {OUTPUT_DIR}: Directory in which to house the RR counts. This will be created if not in existence.
- {NORMED}: Whether to normalize the output RR counts by the summed squared random weights (0 = no normalization, 1 = normalization). This should be done in LEGENDRE mode (where the RR counts are used to compute the survey correction function) but *not* in the DEFAULT mode.

**Notes**:

- Output files will be saved as ``RR_counts_n{N}_m{M}_{INDEX}.txt`` for N radial and M angular bins. INDEX specifies which fields are used (e.g. INDEX=12 for the field 1 x field 2 pair count). This is saved as a simple list of :math:`N\times M` values using the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}*M + \mathrm{bin}_\mathrm{angular}`.
- This is a very CPU intensive computation since we must compute pair counts between every pair of random particles up to the maximum particle separation. The process can be expedited using multiple CPU cores or a reduced number of random particles (e.g. via the :ref:`particle-subset` script).
- For two sets of input particles, three sets of weights must be computed for the three possible pairs of two distinct fields, hence the computation time increases by a factor of three.


.. _RRR_counts:

Estimating RRR triple counts
-----------------------------

.. todo:: write this


.. _survey_correction_2PCF:

Computing Survey Correction Functions
--------------------------------------

This script computes the survey correction function :math:`\Phi(r_a,\mu)` in each bin :math:`a` and fits it to a smooth model. The output parameters can be fed into the main C++ script. For convenience, we divide the output :math:`\Phi` by :math:`V\overline{(nw)^2}`, for survey volume :math:`V`.

Usage
~~~~~~

For a single field analysis::

    python python/compute_correction_function.py {GALAXY_FILE} {BIN_FILE} {RR_counts} {OUTPUT_DIR}
    
For an analysis using two distinct fields::

    python python/compute_correction_function_multi.py {GALAXY_FILE} {GALAXY_FILE_2} {BIN_FILE} {RR_COUNTS_11} {RR_COUNTS_12} {RR_COUNTS_22} {OUTPUT_DIR}

**Input Parameters**:

See the :ref:`RR_counts` parameters above. Additionally;

- {RR_COUNTS}, {RR_COUNTS_11}, {RR_COUNTS_12} {RR_COUNTS_22}: RR counts computed by the :ref:`RR_counts` script. These should be normalized by the sum of the squared random weights (by using the NORMED flag above).

**Notes**:

- Output
.. todo:: test this script and note that we add nw2 V normalizations

.. todo:: add 3PCF to here
