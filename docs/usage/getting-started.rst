Getting Started
================

RascalC computes covariance matrix estimates from a given correlation function and set of random particles. Here, we give a broad overview of the procedure and the relevant :ref:`file-inputs`. A :doc:`tutorial` is provided to demonstrate basic use of the code pipeline.

In order to compute the covariance matrices there are several steps:

1. :doc:`pre-processing` (*Optional*):
    We provide a suite of utility functions to convert input files to the correct forms used by RascalC. This includes conversion from (Ra,Dec,redshift) to (x,y,z) coordinate space, creation of binfiles and assignment of HealPix jackknife regions to particles. Alternatively, this step can be skipped if the input files are already of the correct format.
2. :doc:`jackknife-weights`: 
    Before the main C++ code is run, we compute the weights for each jackknife region, by computing jackknife-specific RR pair counts using `Corrfunc <https://corrfunc.readthedocs.io>`_. This is run via a Python script.
3. :doc:`correlation-functions` (*Optional*): 
    This provides functions to compute the jackknife correlation functions :math:`\xi^{J}(r,\mu)` for one or two input fields (using Corrfunc), which are later used to calibrate the shot-noise rescaling parameter(s). In addition, we provide routines to compute the overall survey correlation functions. These may also be defined by the user instead.
4. :doc:`main-code`:
    The main C++ code computing the individual covariance matrix terms using Monte Carlo integration. For multiple input correlation functions, this computes all relevant terms for the six non-trivial cross-covariance matrices. The covariances are saved as ``.txt`` files which can be reconstructed in Python.
5. :doc:`post-processing`: 
    This Python script computes the shot-noise rescaling parameter(s) and reconstructs output covariance matrices from the jackknive correlation function estimates produced in :doc:`correlation-functions`. A single ``.npz`` file is created including the output covariance and precision matrices as well as the effective number of mocks :math:`N_\mathrm{eff}`.
6. :doc:`visualization`: (*Optional*)
    We provide a Python module for reconstruction and visualization of the output covariance matrices. This performs similar routines to the :doc:`post-processing` codes, but allows for greater user interactivity, as well as plotting regimes.

.. todo:: add in post-processing multi-field scripts.
    
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
    - *Format*: See above.
- *(Usually created internally)* **Correlation Function(s)**:
    - This specifies the input correlation function estimates to be used by RascalC. 
    - For two sets of tracer particles, we require three correlation functions; two auto-correlations and a cross-correlation.
    - These can be user input or created with Corrfunc using the :ref:`full-correlations` codes.
    - Estimates of :math:`\xi(r,\mu)` must be given for a grid of values of :math:`(r,\mu)`, which must extend close to zero for :math:`r` with the bins as specified in the correlation function binning file.
    - *Format*: An ASCII file with space separated values. Line 1 lists the radial coordinates of the bin centers and line 2 lists the angular coordinates. Successive lines list the correlation function estimates :math:`\xi(r,\mu)`, with the column indicating the :math:`\mu` bin center and the row indicating the :math:`r` bin center.
- *(Usually created internally)* **Jackknife Correlation Functions**:
    - This specifies the input correlation function estimates for each *unrestricted* jackknife, :math:`\xi^J_{A}(r,\mu)`. 
    - For two sets of tracer particles, we require three correlation functions; two auto-correlations and a cross-correlation.
    - This is conventionally created with Corrfunc using the :ref:`jackknife-correlations` codes, but may be user input if desired.
    - The radial and angular binning should match that desired for the output covariance matrix.
    - *Format*: An ASCII file with space separated values. Lines 1 and 2 list the radial and angular bin centers (as for the full correlation function). Each succeeding line gives the entire correlation function estimate for a given jackknife. The rows indicate the jackknife and the columns specify the collapsed bin, using the indexing :math:`\mathrm{bin}_\mathrm{collapsed} = \mathrm{bin}_\mathrm{radial}\times n_\mu + \mathrm{bin}_\mathrm{angular}` for a total of :math:`n_\mu` angular bins (unlike for the full correlation function). 
- *(Usually created internally)* **Jackknife Weights and Random Particle Counts**:
    - These specify the weights of each jackknife region for each bin and the random particle counts both for each jackknife, and for the entire survey. 
    - These should be created using the :doc:`jackknife-weights` script.
    - They are saved in ``.dat`` files with the name ``jackknife_weights_n{N}_m{M}_j{J}_{INDEX}.dat``, ``jackknife_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}_{INDEX}.dat`` where N and M specify the number of radial and angular bins respectively and J gives the number of non-empty jackknife regions. INDEX specifies which fields are being used (e.g. 12 specifies the cross-weights between fields 1 and 2).
    
