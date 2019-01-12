Getting Started
================

RascalC computes covariance matrix estimates from a given correlation function and set of random particles. The required input files are described :ref:`below <file-inputs>`.

In order to compute these matrices there are 4 steps:

1. :doc:`pre-processing` (*Optional*):
    We provide a suite of utility functions to convert input files to the correct forms used by RascalC. This includes conversion from (Ra,Dec,redshift) to (x,y,z) coordinate space, creation of binfiles and assignment of HealPix jackknife regions to particles. Alternatively, this step can be skipped if the input files are already of the correct format.
2. :doc:`jackknife-weights`: 
    Before the main C++ code is run, we compute the weights for each jackknife region, by computing jackknife-specific RR pair counts using `Corrfunc <https://corrfunc.readthedocs.io>`_. This is run via a Python script.
3. :doc:`main-code`:
    The main C++ code computing the individual covariance matrix terms using Monte Carlo integration. For multiple input correlation functions, this computes all relevant terms for the six non-trivial cross-covariance matrices. The covariances are saved as ``.txt`` files which can be reconstructed in Python.
4. :doc:`post-processing`: (*Optional*)
    A suite of codes in Python are provided to reconstruct and save the output covariance matrices. In addition, we provide modules to compute precision matrix estimates and effective mock numbers as well as plotting routines.

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
- **Covariance Matrix Binning File**:
    - This specifies the radial binning in the output covariance matrix.
    - For each bin we specify the minimum and maximum radii in comoving Mpc/h units.
    - Linear, logarithmic and hybrid binning files can be created using the :ref:`write-binning-file` scripts.
    - *Format*: An ASCII file with each bin occupying a separate line, with tab-separated columns specifying :math:`(r_\mathrm{min},r_\mathrm{max})` for each bin.
- **Correlation Function(s)**:
    - This specifies the input correlation function estimates to be used by RascalC. 
    - For two sets of tracer particles, we require three correlation functions; two auto-correlations and a cross-correlation.
    - *Coming Soon*: Code to estimate the correlation functions from galaxy positions using Corrfunc.
    - Estimates of :math:`\xi(r,\mu)` must be given for a grid of values of :math:`(r,\mu)`, which must extend close to zero for :math:`r` and cover the entire range of radial bins specified in the binning file.
    - *Format*: An ASCII file with space separated values. Line 1 lists the radial coordinates of the bin centers and line 2 lists the angular coordinates. Successive lines list the correlation function estimates :math:`\xi(r,\mu)`, with the column indicating the :math:`\mu` bin and the row indicating the :math:`r` bin.
    
.. todo:: add correlation function creation script

- *(Internally Created)* **Jackknife Weights and Random Particle Counts**:
    - These specify the weights of each jackknife region for each bin and the random particle counts for each jackknife. 
    - These must be created using the :doc:`jackknife-weights` script.
    - They are saved in ``.dat`` files with the name ``jackknife_weights_n{N}_m{M}_j{J}.dat`` and ``binned_pair_counts_n{N}_m{M}_j{J}.dat`` where N and M specify the number of radial and angular bins respectively and J gives the number of non-empty jackknife regions.
    
.. todo:: add support for multi-tracer jackknife weights
