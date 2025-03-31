Getting Started
================

RascalC computes covariance matrix estimates from a given correlation function and set of random particles. Here, we give a broad overview of the procedure and the relevant :ref:`file-inputs`. To demonstrate basic use of the code pipeline, we provide two tutorials: :doc:`tutorial` and :doc:`tutorial_periodic`.

.. _computation_modes:

Computation Modes
------------------

The RascalC code may be run a variety of modes, which compute the covariances of different statistics. These are enabled via compilation flags, set in the Makefile (see :doc:`main-code` for usage). The various modes are outlined below:

1. **DEFAULT** (No compiler flags): Compute the full-survey covariance of the anisotropic two-point correlation function (2PCF) in :math:`(r,\mu)` bins.
2. **LEGENDRE** (``-DLEGENDRE`` flag): Compute the full-survey covariance of (even) Legendre multipoles of the 2PCF accumulated directly. *NB*: We do not provide functionality to compute the jackknife covariance in Legendre multipole bins, since this has a more complex functional form. However, the shot-noise rescaling parameter :math:`\alpha` can be found from :math:`(r,\mu)` jackknife covariance matrix fitting and applied to the output full-survey Legendre-binned matrix
3. **LEGENDRE_MIX** (``-DLEGENDRE_MIX`` flag): Compute the full-survey covariance of (even) Legendre multipoles of the 2PCF projected from :math:`(r,\mu)` bins. Compatible with jackknives.
4. **JACKKNIFE** (``-DJACKKNIFE`` flag): Compute the full-survey and jackknife covariance matrices of the anisotropic 2PCF in :math:`(r,\mu)` bins. The theoretical jackknife matrix can be compared to the sample jackknife matrix to compute the shot-noise rescaling parameter :math:`\alpha`. Compatible with **DEFAULT** and **LEGENDRE_MIX** modes but neither **LEGENDRE** nor **3PCF**.
5. **3PCF** (``-DTHREE_PCF`` flag): Compute the full-survey covariance of (odd and even) Legendre multipoles of the isotropic three-point correlation function (3PCF).

.. _pipeline_outline:

Outline of Pipeline
--------------------

In order to compute the covariance matrices there are several steps:

1. :doc:`pre-processing` (*Optional*):
    We provide a suite of utility functions to convert input files to the correct forms used by RascalC. This includes conversion from (Ra,Dec,redshift) to (x,y,z) coordinate space, creation of binfiles and assignment of HealPix jackknife regions to particles. Alternatively, this step can be skipped if the input files are already of the correct format.
2. :doc:`jackknife-weights` (*Only required in JACKKNIFE mode*):
    Before the main C++ code is run, we compute the weights for each jackknife region, by computing jackknife-specific RR pair counts using `Corrfunc <https://corrfunc.readthedocs.io>`_. This is run via a Python script.
3. 
    a. :doc:`geometry-correction` (*Required in DEFAULT, LEGENDRE, LEGENDRE_MIX and 3PCF modes*):
        The main C++ code requires either the RR counts (DEFAULT mode) or the survey-geometry correction function :math:`\Phi` (LEGENDRE and 3PCF modes). We provide Python and C++ scripts to compute these, with the correction functions requiring a set of RR or RRR counts to be computed first.
    b. :doc:`correlation-functions` (*Optional*):
        This provides functions to compute the overall survey correlation functions for one or two fields using Corrfunc (which may also be defined by the user). In addition, we provide routines to compute the jackknife correlation functions :math:`\xi^{J}(r,\mu)`, which are later used to calibrate the shot-noise rescaling parameter(s) in the JACKKNIFE mode.
4. :doc:`main-code`:
    The main C++ code computing the 2PCF or 3PCF individual covariance matrix terms using Monte Carlo integration. For multiple input correlation functions in the 2PCF modes, this computes all relevant terms for the six non-trivial cross-covariance matrices. The covariances are saved as ``.txt`` files which can be reconstructed in Python.
5. :doc:`post-processing`:
    If jackknife covariances have been computed, this Python script computes the shot-noise rescaling parameter(s) and reconstructs output covariance matrices from the jackknive correlation function estimates produced in :doc:`correlation-functions`. Else, we provide scripts to reconstruct the :doc:`main-code` output. In both cases, a single ``.npz`` file is created including the output covariance and precision matrices as well as the effective number of mocks :math:`N_\mathrm{eff}`.

.. _file-inputs:

File Inputs
------------

The required input files and formats are described below. Note that several of these can be computed using the :doc:`pre-processing` codes.

- **Random Particle File(s)**:
    - This lists the locations and weights of random particles which describe a survey geometry.
    - This must specify the {x,y,z,w,j} coordinates for each particle, where {x,y,z} are Cartesian coordinates (in comoving Mpc/h units), w are particle weights and j are integers referencing which jackknife the particle is in.
    - {RA,Dec,redshift} coordinates can be converted to {x,y,z} positions using the :ref:`coord-conversion` script.
    - HealPix jackknives can be added using the :ref:`create-jackknives` script.
    - *Format*: An ASCII file with each particle defined on a new row, and tab-separated columns indicating the {x,y,z,w,j} coordinates.
- **Galaxy Position File(s)**:
    - This lists the locations and weights of galaxies in a specific survey, in the same manner as the random particles.
    - This is only required to compute the correlation functions in the :doc:`correlation-functions` scripts.
    - *Format*: See above.
- **Covariance Matrix Binning File**:
    - This specifies the radial binning in the output covariance matrix.
    - For each bin we specify the minimum and maximum radii in comoving Mpc/h units.
    - Linear, logarithmic and hybrid binning files can be created using the :ref:`write-binning-file` scripts.
    - *Format*: An ASCII file with each bin occupying a separate line, with tab-separated columns specifying :math:`(r_\mathrm{min},r_\mathrm{max})` for each bin.
- **Correlation Function Binning File**:
    - File specifying the radial binning used in the input correlation function.
    - The lowest bin must extend to zero for this, and the highest bin should be at least as large as the maximum covariance matrix bin.
    - Currently must be the same for all input correlation functions, for the multiple field case.
    - *Format*: See above.
- *(Usually created internally)* **Correlation Function(s)**:
    - This specifies the input correlation function estimates to be used by RascalC.
    - For two sets of tracer particles, we require three correlation functions; two auto-correlations and a cross-correlation.
    - These can be user input or created with Corrfunc using the :ref:`full-correlations` codes.
    - Estimates of :math:`\xi(r,\mu)` must be given for a grid of values of :math:`(r,\mu)`, which must extend close to zero for :math:`r` with the bins as specified in the correlation function binning file.
    - *Format*: An ASCII file with space separated values. Line 1 lists the radial coordinates of the bin centers and line 2 lists the angular coordinates. Successive lines list the correlation function estimates :math:`\xi(r,\mu)`, with the column indicating the :math:`\mu` bin center and the row indicating the :math:`r` bin center.
- *(Only required in JACKKNIFE mode and usually created internally)* **Jackknife Correlation Functions**:
    - This specifies the input correlation function estimates for each *unrestricted* jackknife, :math:`\xi^J_{A}(r,\mu)`.
    - For two sets of tracer particles, we require three correlation functions; two auto-correlations and a cross-correlation.
    - This is conventionally created with Corrfunc using the :ref:`jackknife-correlations` codes, but may be user input if desired.
    - The radial and angular binning should match that desired for the output covariance matrix.
    - If this is supplied separately, the user must ensure that the pair count terms are normalized by the ratio of summed galaxy and random particle weights across the **entire** survey, not just those in the relevant jackknife region. This is for later convenience when estimating the jackknife covariance matrix model.
    - *Format*: An ASCII file with space separated values. Lines 1 and 2 list the radial and angular bin centers (as for the full correlation function). Each succeeding line gives the entire correlation function estimate for a given jackknife. The rows indicate the jackknife and the columns specify the collapsed bin, using the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins (unlike for the full correlation function).
- *(Required in JACKKNIFE mode with DEFAULT or LEGENDRE_MIX and usually created internally)* **Jackknife Weights and Random Particle Counts**:
    - These specify the weights of each jackknife region for each bin and the random particle counts both for each jackknife, and for the entire survey.
    - These should be created using the :doc:`jackknife-weights` script.
    - They are saved in ``.dat`` files with the name ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat``, ``jackknife_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` where N and M specify the number of radial and angular bins respectively and J gives the number of non-empty jackknife regions. INDEX specifies which fields are being used (e.g. 12 specifies the cross-weights between fields 1 and 2).
- *(Required in DEFAULT or LEGENDRE_MIX mode and usually created internally)* **Random Particle Counts**:
    - These specify random particle counts for the entire survey, which are needed to normalize the :math:`(r,\mu)` binned covariances.
    - These should be created using the RR count script described in :doc:`geometry-correction` (and *not* normalized by the summed squared weights).
    - They are saved in ``.dat`` files with the name ``binned_pair_counts_n{N}_m{M}_{INDEX}.dat`` where N and M specify the number of radial and angular bins respectively. INDEX specifies which fields are being used (e.g. 12 specifies the cross-weights between fields 1 and 2).
- *(Required in LEGENDRE and 3PCF modes and usually created internally)* **Survey Correction Function Parameters**:
    - These give the necessary parameters for the main C++ code to reconstruct the survey-correction function, :math:`\Phi(r_a,\mu)` (2PCF) or :math:`\Phi(r_a,r_b,\chi)` (3PCF).
    - For multiple input fields, we will have three output bin correction factors of the same format.
    - These should be created using the survey-correction functions described in :doc:`geometry-correction`, and require the RR or RRR counts to be computed (also described in :doc:`geometry-correction`).
    - They are saved as ASCII files with the names ``BinCorrectionFactor_n{N}_m{M}.txt`` or ``BinCorrectionFactor3PCF_n{N}_m{M}.txt`` and specify polynomial fitting parameters (2PCF) or the first seven multipoles of :math:`\Phi^{-1}` (3PCF), which are found to well describe the fit. These have one row per radial bin (or pair of bins for the 3PCF), and must be constructed using the same radial binning as for the output covariance matrix.
- *(Required in LEGENDRE_MIX mode and usually created internally)* **Projection factors from** :math:`\mu` **bins to Legendre multipoles**:
    - These together with the full-survey random particle counts give the correct normalization for the projected Legendre multipole covariance.
    - One file of them is enough, it can be created with the :ref:`mu_bin_legendre_factors` script.
    - The file must have rows corresponding to the :math:`\mu` bins and columns corresponding to the (even) Legendre multipoles. The factors are the same for all radial bins, unlike the random counts which also influence the projection.
