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
    
*(See :ref:`coord-conversion`)*.
    
These are saved as ``qpm_galaxies.xyz`` and ``qpm_randoms_50x.xyz`` in (x,y,z,weight) format.

Now let's add some jackknives to these files. We'll use HEALPIX NSIDE=8 jackknife regions here::

    python python/create_jackknives.py qpm_galaxies.xyz qpm_galaxies.xyzwj 8
    python python/create_jackknives.py qpm_randoms_50x.xyz qpm_randoms_50x.xyzwj 8

*(See :ref:`create-jackknives`)*.

These are saved as ``qpm_galaxies.xyzwj`` and ``qpm_randoms_50x.xyzwj`` in (x,y,z,weight,jackknife-id) format, and we're using 169 non-empty jackknives here.

We've got 50x the number of random particles as galaxies here which seems a little big. Let's reduce this to 10x (noting that there are 642051 galaxies in the galaxy file)::

    python python/take_subset_of_particles.py qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj 6420510
    
*(See :ref:`particle-subset`)*.
    
Great! Now we have a random particle file with 10x randoms, and all files are in the correct format. 

Let's create the radial binning files next. We'll create two binning files; one for the correlation functions and one for the covariance matrix.

For the covariance matrix, we'll use a linear binning file with :math:`\Delta r = 5` for :math:`r\in[20,200]` and for the correlation function we'll use a linear binning file with :math:`\Delta r = 1` for :math:`r\in[0,200]`. **NB**: The correlation function binning file must extend down to :math:`r = 0`::

    python python/write_binning_file_linear.py 36 20 200 radial_binning_cov.csv
    python python/write_binning_file_linear.py 200 0 200 radial_binning_corr.csv
    
*(See :ref:`write-binning-file`)*.

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

*(See :doc:`jackknife-weights`)*.

This computes pair counts for each pair of random particles in the survey (using Corrfunc), so may take a while...

The outputs are saved as ``weights/jackknife_weights_n36_m12_j169_11.dat``, ``weights/jackknife_pair_counts_n36_m12_j169_11.dat`` and ``weights/binned_pair_counts_n36_m12_j169_11.dat`` containing the weights :math:`w_{aA}`, bin-counts :math:`RR_{aA}` and summed bin counts :math:`RR_a` respectively.

3) Correlation Functions
-------------------------

Using the galaxy and random particle files, we can obtain estimates of the correlation function. Firstly, we'll compute an estimate of :math:`\xi(r,\mu)` to be used to compute the theoretical covariance matrices. 
In addition, we'll use 120 :math:`\mu` bins in :math:`[0,1]` and set the code to run for aperiodic input data. This must use the *correlation function* radial binning file, giving us a fine estimate of the correlation function.::

    python python/xi_estimator.py qpm_galaxies.xyzwj qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj radial_binning_corr.csv 1. 120 10 0 xi/
    
*(See :ref:`full-correlations`)*

This uses Corrfunc to perform pair counting and computes :math:`\xi(r,\mu)` via the Landy-Szalay estimator. Here we're using 10x randoms to compute the RR pair counts and 50x randoms to compute the DR pair counts. The output is saved as ``xi/xi_n200_m120_11.dat`` in the format specified in :ref:`file-inputs`.



Do stuff::

    python python/xi_estimator_jack.py {GALAXY_FILE_1} {GALAXY_FILE_2} {RANDOM_FILE_1_DR} {RANDOM_FILE_1_RR} {RANDOM_FILE_2_DR} {RANDOM_FILE_2_RR} {RADIAL_BIN_FILE} {MU_MAX} {N_MU_BINS} {NTHREADS} {PERIODIC} {OUTPUT_DIR} [{RR_counts_11} {RR_counts_12} {RR_counts_22}]

*(See :ref:`jackknife-correlations`)*

