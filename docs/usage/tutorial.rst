Usage Tutorial
===============


We present a basic example of the use of the RascalC code for a single field. Multiple field cases proceed similarly. Detailed documentation for all functions is given in associated pages, as overviewed in the :doc:`getting-started` pages.


Here, we compute the covariance matrix for a single `QPM <https://arxiv.org/pdf/1309.5532.pdf>`_ mock dataset. We'll work in the directory in which RascalC is installed for simplicity.

1) Pre-Processing
------------------

**Inputs**:

- ``mock_galaxy_DR12_CMASS_N_QPM_0001.txt``: Galaxy positions and weights for the QPM mock in (Ra,Dec,redshift,weight) format.
- ``mock_random_DR12_CMASS_N_50x1.txt``: 50x Random positions and weights for the CMASS survey in (Ra,Dec,redshift,weight).

First, we'll convert these into Cartesian (x,y,z,weight) coordinates, using :math:`\Omega_m = 0.31`, :math:`\Omega_k = 0`, :math:`w_\Lambda = -1` (to be consistent with BOSS DR12)::

    python python/convert_to_xyz.py mock_galaxy_DR12_CMASS_N_QPM_0001.txt qpm_galaxies.xyz 0.31 0. -1
    python python/convert_to_xyz.py mock_random_DR12_CMASS_N_50x1.txt qpm_randoms_50x.xyz 0.31 0. -1
    
(See :ref:`coord-conversion`).
    
These are saved as ``qpm_galaxies.xyz`` and ``qpm_randoms_50x.xyz`` in (x,y,z,weight) format.

Now let's add some jackknives to these files. We'll use HEALPIX NSIDE=8 jackknife regions here::

    python python/create_jackknives.py qpm_galaxies.xyz qpm_galaxies.xyzwj 8
    python python/create_jackknives.py qpm_randoms_50x.xyz qpm_randoms_50x.xyzwj 8

(See :ref:`create-jackknives`).

These are saved as ``qpm_galaxies.xyzwj`` and ``qpm_randoms_50x.xyzwj`` in (x,y,z,weight,jackknife-id) format, and we're using 169 non-empty jackknives here.

We've got 50x the number of random particles as galaxies here which seems a little big. Let's reduce this to 10x (noting that there are 642051 galaxies in the galaxy file)::

    python python/take_subset_of_particles.py qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj 6420510
    
(See :ref:`particle-subset`).
    
Great! Now we have a random particle file with 10x randoms, and all files are in the correct format. 

Let's create the radial binning files next. We'll create two binning files; one for the correlation functions and one for the covariance matrix.

For the covariance matrix, we'll use a linear binning file with :math:`\Delta r = 5` for :math:`r\in[20,200]` and for the correlation function we'll use a linear binning file with :math:`\Delta r = 1` for :math:`r\in[0,200]`. **NB**: The correlation function binning file must extend down to :math:`r = 0`::

    python python/write_binning_file_linear.py 36 20 200 radial_binning_cov.csv
    python python/write_binning_file_linear.py 200 0 200 radial_binning_corr.csv
    
(See :ref:`write-binning-file`).

Here we're using 36 radial bins for the covariance matrix. Let's have a look at the covariance binning file::

    20.00000000     25.00000000
    25.00000000     30.00000000
    30.00000000     35.00000000
    35.00000000     40.00000000
    40.00000000     45.00000000
    45.00000000     50.00000000
    50.00000000     55.00000000
    55.00000000     60.00000000
    ....
    
This all looks as expected.


2) Jackknife Weights
----------------------

We're now ready to compute the jackknife weights :math:`w_{aA}` for this set of random particles. This determines how much weight we assign to each jackknife region later in the analysis, via the :math:`RR` pair counts in each bin and jackknife.

Here we'll use 12 angular bins with :math:`\mu\in[0,1]` and recall that this dataset is non-periodic (so :math:`\mu` is measured from the line-of-sight, as opposed to the :math:`z`-axis). We'll run 10-threaded for speed and save in the ``weights/`` directory.::

    python python/jackknife_weights.py qpm_randoms_10x.xyzwj radial_binning_cov.csv 1. 12 10 0 weights/

(See :doc:`jackknife-weights`).

This computes pair counts for each pair of random particles in the survey (using Corrfunc), so may take a while...

The outputs are saved as ``weights/jackknife_weights_n36_m12_j169_11.dat``, ``weights/jackknife_pair_counts_n36_m12_j169_11.dat`` and ``weights/binned_pair_counts_n36_m12_j169_11.dat`` containing the weights :math:`w_{aA}`, bin-counts :math:`RR_{aA}` and summed bin counts :math:`RR_a` respectively.


3) Correlation Functions
-------------------------

Using the galaxy and random particle files, we can obtain estimates of the correlation function. Firstly, we'll compute an estimate of :math:`\xi(r,\mu)` to be used to compute the theoretical covariance matrices. 
In addition, we'll use 120 :math:`\mu` bins in :math:`[0,1]` and set the code to run for aperiodic input data. This must use the *correlation function* radial binning file, giving us a fine estimate of the correlation function.::

    python python/xi_estimator.py qpm_galaxies.xyzwj qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj radial_binning_corr.csv 1. 120 10 0 xi/
    
(See :ref:`full-correlations`).

This uses Corrfunc to perform pair counting and computes :math:`\xi_a` for each bin, :math:`a`, via the Landy-Szalay estimator. Here we're using 10x randoms to compute the RR pair counts and 50x randoms to compute the DR pair counts. The output is saved as ``xi/xi_n200_m120_11.dat`` in the format specified in :ref:`file-inputs`. We'll use this full correlation function to compute the theoretical covariance matrix later on.

Now let's compute the jackknnife correlation function estimates for each bin, :math:`\xi^J_{aA}`. These are the individual correlation functions obtained from each unrestricted jackknife, and we can use them to create a data jackknife covariance matrix which we can compare to theory. This is run in a similar way to before, but we must now use the *covariance matrix* radial binning file, since we use these to directly compute a covariance. Here, we'll use 10x randoms for RR counts and 50x randoms for DR counts, but we can skip some of the work by loading in the jackknife pair counts computed by the :doc:`jackknife-weights` script (in the same binning as here), which avoids recomputing RR counts. (The input 10x random file isn't loaded in this case).::

    python python/xi_estimator_jack.py qpm_galaxies.xyzwj qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj radial_binning_cov.csv 1. 12 10 0 xi_jack/ weights/jackknife_pair_counts_n36_m12_j169_11.dat

(See :ref:`jackknife-correlations`).

**NB**: This may take a little while to compute, depending on the number of randoms and galaxies used. The output jackknife correlation functions are saved as ``xi_jack/xi_jack_n36_m12_j169_11.dat`` in the format specified in :ref:`file-inputs`. These will be automatically read later on.


4) Computing the Covariance Matrix
------------------------------------

(See :doc:`main-code`).

Now that all of the inputs have been computed, we can run the main C++ code to compute the theoretical covariance matrix terms. 

There's two ways to run the code here; firstly we could edit parameters in the ``modules/parameters.h`` file, to tell the code where to find the relevant inputs. Here are the important lines::

    ....
    
    //---------- ESSENTIAL PARAMETERS -----------------
    
    // The name of the input random particle files (first set)
    char *fname = NULL;
    const char default_fname[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/qpm_randoms_10x.xyzwj"; 
    
    // Name of the radial binning .csv file
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/radial_binning_cov.csv";
    
    // The name of the correlation function file for the first set of particles
    char *corname = NULL;
    const char default_corname[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/xi/xi_n200_m120_11.dat";
    
    // Name of the correlation function radial binning .csv file
    char *radial_bin_file_cf = NULL;
    const char default_radial_bin_file_cf[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/radial_binning_corr.csv";
    
    // Number of galaxies in first dataset
    Float nofznorm=642051;
    
    // Name of the jackknife weight file
    char *jk_weight_file = NULL; // w_{aA}^{11} weights
    const char default_jk_weight_file[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/weights/jackknife_weights_n36_m12_j169_11.dat";
    
    // Name of the RR bin file
    char *RR_bin_file = NULL; // RR_{aA}^{11} file
    const char default_RR_bin_file[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/weights/binned_pair_counts_n36_m12_j169_11.dat";
    
    // Output directory 
    char *out_file = NULL;
    const char default_out_file[500] = "/mnt/store1/oliverphilcox/Mock1QPM2/";
    
    // The number of mu bins
    int mbin = 12;
    
    // The number of mu bins in the correlation function
    int mbin_cf = 120;
    
    // The number of threads to run on
    int nthread=10;

    // The grid size, which should be tuned to match boxsize and rmax. 
    // This uses the maximum width of the cuboidal box.
    int nside = 251;

    ....
    
    //---------- PRECISION PARAMETERS ---------------------------------------
        
    // Maximum number of iterations to compute the C_ab integrals over
    int max_loops=10;
    
    // Number of random cells to draw at each stage
    int N2 = 20; // number of j cells per i cell
    int N3 = 40; // number of k cells per j cell
    int N4 = 80; // number of l cells per k cell
    
    ....
    
Here we're using 10 loops (to get 10 independent estimates of the covariance matrix), and setting N2-N4 such that we'll get good precision in a few hours of runtime. Now, we'll compile the code;::
    
    bash clean
    make
 
The first line simply cleans the pre-existing ``./cov`` file, if present and the second compiles ``grid_covariance.cpp`` using the Makefile (using the g++ compiler by default). If we were using periodic data we'd need to set the ``-DPERIODIC`` flag in the Makefile before running this step. Similarly, we could remove the ``-DOPENMP`` flag to run single threaded. The code is then run with the default parameters;

.. code-block:: bash

    ./cov -def
    
Alternatively, we could simply pass these arguments on the command line (after the code is compiled). (**NB**: We can get a summary of the inputs by simply running ``./cov`` with no parameters)

.. code-block:: bash

    ./cov -in qpm_randoms_10x.xyzwj -binfile radial_binning_cov.csv -cor xi/xi_n200_m120_11.dat -binfile_cf radial_binning_corr.csv -norm 642051 -jackknife weights/jackknife_pair_counts_n36_m12_j169_11.dat -RRbin weights/binned_pair_counts_n36_m12_j169_11.dat -output ./ -mbin 12 -mbin_cf 120 -nside 251 -maxloops 10 -N2 20 -N3 40 -N4 80
    
It's often just easier to edit the ``modules/parameter.h`` file, but the latter approach allows us to change parameters without recompiling the code.

This runs in around 5 hours on 10 cores here, giving output matrix components saved in the ``CovMatricesFull`` and ``CovMatricesJack`` directories as ``.txt`` files. We'll now reconstruct these.


5) Post-Processing
-------------------

Although the C++ code computes all the relevant parts of the covariance matrices, it doesn't perform any reconstruction, since this is much more easily performed in Python. Post-processing is used to compute the optimal value of the shot-noise rescaling parameter :math:`\alpha` (by comparing the data-derived and theoretical covariance matrices), as well as construct the output covariance and precision matrices.

For a single field analysis, this is run as follows, specifying the jackknife correlation functions, output covariance term directory and weights. Since we used :math:`N_\mathrm{loops}=10` above, we'll set this as the number of subsamples here::

    python python/post_process.py xi_jack/xi_jack_n36_m12_j169_11.dat weights/ ./ 12 10 ./

(See :ref:`post-processing-single`).
    
The output is a single compressed Python ``.npz`` file which contains the following analysis products:
    - Optimal shot-noise rescaling parameter :math:`\alpha^*`
    - Full theory covariance matrix :math:`C_{ab}(\alpha^*)`
    - Jackknife theory covariance matrix :math:`C^J_{ab}(\alpha^*)`
    - Jackknife data covariance matrix :math:`C^{J,\mathrm{data}}_{ab}
    - Full (quadratic bias corrected) precision matrix :math:`\Psi_{ab}(\alpha^*)`
    - Jackknife (quadratic bias corrected) precision matrix :math:`\Psi^J_{ab}(\alpha^*)`
    - Full quadratic bias :math:`\hat{D}_{ab}` matrix
    - Effective number of mocks :math:`N_\mathrm{eff}`
    - Individual full covariance matrix estimates :math:`C_{ab}^{(i)}(\alpha^*)`
    
This completes the analysis!

.. todo:: link API for visualization??
