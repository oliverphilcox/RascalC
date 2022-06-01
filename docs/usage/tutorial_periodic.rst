Tutorial: Periodic Data and Legendre Multipoles
================================================

We present a basic example of the use of the RascalC code on a periodic dataset, such as the output of an N-body simulation. Here we'll restrict to a single-field analysis (using only a single set of tracer galaxies), but the multiple field case proceeds similarly. Our goal is to compute the covariance matrix of the monopole and quadrupole correlation functions :math:`\xi_0(r)` and :math:`\xi_2(r)`, as described in `Philcox & Eisenstein 2019 <https://arxiv.org/abs/1910.04764>`_. For this we will use the code in LEGENDRE mode. Detailed documentation for all functions is given in associated pages, as overviewed in the :doc:`getting-started` pages. Note that we here run the main code in the JACKKNIFE mode, and work in the RascalC installation directory for simplicity.

An additional tutorial (:doc:`tutorial`) shows the usage of RascalC for an *aperiodic* dataset (e.g. a galaxy survey) in :math:`(r,\mu)` co-ordinates. This also includes descriptions of the shot-noise rescaling procedure used to incorporate non-Gaussianity in our covariances. In this tutorial, we will assume full Gaussianity, setting the shot-noise rescaling parameter to :math:`\alpha=1`. To include non-Gaussianity, we simply need to run RascalC in JACKKNIFE mode first to estimate :math:`\alpha` (as in :doc:`tutorial`), then run it in LEGENDRE mode, using the measured value of :math:`\alpha` in post-processing.

1) Pre-Processing
------------------

**Inputs**:

- ``nbody_simulation.txt``: Galaxy positions and weights for the periodic mock dataset in (x,y,z,weight) format (in Mpc/h). We assume this to have a cubic-box geometry.

First, we will create a set of *random particles* in a box of the same shape as our simulation data-set. These will be used as sampling points for our covariance matrix estimator. The simulation had boxsize :math:`L = 1000\,h^{-1}\mathrm{Mpc}` with :math:`N_\mathrm{gal}=156800` galaxies, so we'll create a random set with the same boxsize and :math:`N_\mathrm{rand} = 10 N_\mathrm{gal}` particles. This can be done with the simple Python script

.. code-block:: python

    import numpy as np
    boxsize = 1000 # boxsize in Mpc/h
    N_rand = 1568000 # number of random particles

    # Define random positions and (uniform) weights
    rand_position = np.random.uniform(low=0,high=1,size=(N_rand,3))*boxsize
    rand_weight = np.ones_like(rand_position[:,0]).reshape(-1,1)

    # Concatenate data:
    rand_data = np.hstack([rand_position,rand_weight])

    # Save to file
    np.savetxt('nbody_randoms_10x.txt',rand_data)

This is particularly easy here since N-body simulations have uniform selection functions and simple geometries.

Let's create the radial binning files next. We'll create two binning files; one for the correlation functions and one for the covariance matrix.

For the covariance matrix, we'll use a linear binning file with :math:`\Delta r = 5` for :math:`r\in[25,150]` and for the correlation function we'll use a linear binning file with :math:`\Delta r = 1` for :math:`r\in[0,200]`. **NB**: The correlation function binning file must extend down to :math:`r = 0`::

    python python/write_binning_file_linear.py 25 25 150 radial_binning_cov.csv
    python python/write_binning_file_linear.py 200 0 200 radial_binning_corr.csv

(See :ref:`write-binning-file`).

Here we're using 25 radial bins for the covariance matrix. Let's have a look at the covariance binning file::

    25.00000000     30.00000000
    30.00000000     35.00000000
    35.00000000     40.00000000
    40.00000000     45.00000000
    45.00000000     50.00000000
    50.00000000     55.00000000
    55.00000000     60.00000000
    ....

This all looks as expected.


2) Correlation Functions and Geometry Correction
-------------------------------------------------

Using the galaxy position file, we can obtain estimates of the 2-point correlation function. This is needed to compute the theoretical covariance matrices. Even though our desired output is a covariance matrix of the 2PCF in Legendre multipoles, the main C++ code uses an input 2PCF in :math:`(r,\mu)` space for speed. Here, let's use 20 :math:`\mu` bins in :math:`[0,1]` and the correlation function binning specified above. Since we have periodic data, we must use the ``xi_estimator_periodic.py`` script to compute this. (This is much faster than for the aperiodic case since we can compute the random particle counts analytically.)::

    python python/xi_estimator_periodic.py nbody_simulation.txt radial_binning_corr.csv 1000 1. 20 10 xi/

(See :ref:`full-correlations`).

This uses Corrfunc to perform pair counting and computes :math:`\xi_a` for each bin, :math:`a`, via the standard :math:`DD/RR-1` estimator. Here we're running on 10 cores, assuming a box-size of :math:`L = 1000\,h^{-1}\mathrm{Mpc}`, and the output is saved as ``xi/xi_n200_m20_periodic_11.dat`` in the format specified in :ref:`file-inputs`. We'll use this full correlation function to compute the theoretical covariance matrix later on.

The main C++ code requires an input *survey correction function* to account for non-trivial survey geometries. For a periodic box, the correction function :math:`\Phi(r_a,\mu)` is constant, but the normalization carries important information including the survey volume and number density. This is simply computed via::

    python python/compute_correction_function.py nbody_simulation.txt radial_binning_cov.csv ./ 1

(See :doc:`geometry-correction`)

The :math:`1` specifies a periodic survey, and the code produces the output file ``BinCorrectionFactor_n25_periodic_11.txt`` in the working directory that can be fed into the C++ code. Note that we use the *covariance matrix* binning file here.

3) Computing the Covariance Matrix
------------------------------------

(See :doc:`main-code`).

Now that the necessary inputs have been computed, we can run the main C++ code to compute the theoretical covariance matrix terms.

There's two ways to run the code here; firstly we could edit parameters in the ``modules/parameters.h`` file, to tell the code where to find the relevant inputs. Here are the important lines

.. code-block:: c++

    ....

    //---------- ESSENTIAL PARAMETERS -----------------

    // The name of the input random particle files (first set)
    char *fname = NULL;
    const char default_fname[500] = "nbody_randoms_10x.txt";

    // Name of the radial binning .csv file
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "radial_binning_cov.csv";

    // The name of the correlation function file for the first set of particles
    char *corname = NULL;
    const char default_corname[500] = "xi/xi_n200_m20_periodic_11.dat";

    // Name of the correlation function radial binning .csv file
    char *radial_bin_file_cf = NULL;
    const char default_radial_bin_file_cf[500] = "radial_binning_corr.csv";

    // Number of galaxies in first dataset
    Float nofznorm=156800;

    // Output directory
    char *out_file = NULL;
    const char default_out_file[500] = "./";

    // The number of mu bins in the correlation function
    int mbin_cf = 20;

    // The number of threads to run on
    int nthread=10;

    // The grid size, which should be tuned to match boxsize and rmax.
    // This uses the maximum width of the cuboidal box.
    int nside = 151;

    // Whether or not we are using a periodic box
    bool perbox = true;

    ....

    //-------- LEGENDRE PARAMETERS -------------------------------------------

    int max_l = 2; // max Legendre moment (must be even unless computing 3PCF)

    char *phi_file = NULL; // Survey correction function coefficient file
    const char default_phi_file[500] = "BinCorrectionFactor_n25_periodic_11.txt";

    ....

    //---------- PRECISION PARAMETERS ---------------------------------------

    // Maximum number of iterations to compute the C_ab integrals over
    int max_loops=10;

    // Number of random cells to draw at each stage
    int N2 = 10; // number of j cells per i cell
    int N3 = 10; // number of k cells per j cell
    int N4 = 10; // number of l cells per k cell

    ....

Here we're using 10 loops (to get 10 independent estimates of the covariance matrix), and setting N2-N4 such that we'll get good precision in a few hours of runtime. Now, we'll compile the code;::

    bash clean
    make

The first line simply cleans the pre-existing ``./cov`` file, if present and the second compiles ``grid_covariance.cpp`` using the Makefile (using the g++ compiler by default). We have edited the Makefile to add the ``-DPERIODIC`` flag and ``-DLEGENDRE`` flags to ensure we compute covariances of 2PCF Legendre moments in a periodic geometry. Note that we can also remove the ``-DOPENMP`` flag to run single threaded. The code is then run with the default parameters;

.. code-block:: bash

    ./cov -def

Alternatively, we could simply pass these arguments on the command line (after the code is compiled). (**NB**: We can get a summary of the inputs by simply running ``./cov`` with no parameters)

.. code-block:: bash

    ./cov -in nbody_randoms_10x.txt -binfile radial_binning_cov.csv -cor xi/xi_n200_m20_periodic_11.dat -binfile_cf radial_binning_corr.csv -norm 156800 -output ./ -mbin_cf 20 -nthread 10 -perbox -max_l 2 -phi_file BinCorrectionFactor_n25_periodic_11.txt -maxloops 10 -N2 10 -N3 20 -N4 40 -cf_loops 0

It's often just easier to edit the ``modules/parameter.h`` file, but the latter approach allows us to change parameters without recompiling the code. *NB*: Sometimes, this will crash with the error `Average particle density exceeds maximum advised particle density'; this is due to the ``-nside`` parameter being too low. To fix this, increase ``-nside`` to a larger (odd) value (default; 71).

This runs in under an hour on 10 cores here, giving output matrix components saved in the ``CovMatricesFull`` directory as ``.txt`` files. We'll now reconstruct these.

4) Post-Processing
-------------------

Although the C++ code computes all the relevant parts of the covariance matrices, it doesn't perform any reconstruction, since this is much more easily performed in Python. Post-processing is used to add together the relevant matrix outputs, constructing the full covariance and precision matrices.

For a single field analysis, this is run as follows, specifying the jackknife correlation functions, output covariance term directory and weights::

    python python/post_process_legendre.py ./ 25 2 10 ./ 1.

(See :ref:`post-processing-general`).

Here the first parameter gives the directory where the RascalC products are stored (i.e. the location of the ```CovMatricesFull`` directory), whilst the remained specify number of radial bins, maximum Legendre multipole, number of matrix estimates computed (equal to the ``-maxloops`` parameter in the C++ code) and the output directory. The *optional* last parameter is the shot-noise rescaling parameter. This cannot be computed from Legendre multipole binned data; to include it we must run the code in JACKKNIFE mode first then use the derived shot-noise rescaling parameter in the LEGENDRE mode post-processing (see :doc:`tutorial` for a tutorial on JACKKNIFE mode computations.)

If the covariance matrix terms have not converged sufficiently well (leading to a non-invertable covariance matrix) the Python code will exit prematurely and *not* give a reconstructed output file. This indicates that the main C++ code should be run for longer or with a larger number of random particles.

The output is a single compressed Python ``.npz`` file named ``Rescaled_Covariance_Matrices_Legendre_n25_l2.npz`` which contains the following analysis products:

    - Full theory covariance matrix :math:`C_{ab}(\alpha^*)`
    - Utilized shot-noise rescaling parameter :math:`\alpha^*` (either user-input or set to unity)
    - Effective number of mocks :math:`N_\mathrm{eff}`
    - Full quadratic bias :math:`\tilde{D}_{ab}` matrix
    - Individual full covariance matrix estimates :math:`C_{ab}^{(i)}(\alpha^*)`

Each matrix is stored as a two-dimensional array, with each index specifying both the radial bin and Legendre multipole. The elements are ordered first by Legendre multipole then by radial bin, such that the first elements is the :math:`\ell=0` element in the first radial bin, the second is the :math:`\ell=2` element in the first radial bin, etc.

This completes the analysis!
