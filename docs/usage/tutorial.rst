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
    
(See :ref:`coord-conversion` for documentation).
    
These are saved as ``qpm_galaxies.xyz`` and ``qpm_randoms_50x.xyz`` in (x,y,z,weight) format.

Now let's add some jackknives to these files. We'll use HEALPIX NSIDE=8 jackknife regions here::

    python python/create_jackknives.py qpm_galaxies.xyz qpm_galaxies.xyzwj 8
    python python/create_jackknives.py qpm_randoms_50x.xyz qpm_randoms_50x.xyzwj 8

(See :ref:`create-jackknives` for documentation).

These are saved as ``qpm_galaxies.xyzwj`` and ``qpm_randoms_50x.xyzwj`` in (x,y,z,weight,jackknife-id) format.

We've got 50x the number of random particles as galaxies here which seems a little big. Let's reduce this to 10x (noting that there are 642051 galaxies in the galaxy file)::

    python python/take_subset_of_particles.py qpm_randoms_50x.xyzwj qpm_randoms_10x.xyzwj 6420510
    
(See :ref:`particle-subset` for documentation).
    
Great! Now we have a random particle file with 10x randoms, and all files are in the correct format. 

Let's create the radial binning files next. We'll create two binning files; one for the correlation functions and one for the covariance matrix.

For the covariance matrix, we'll use a linear binning file with :math:`\Delta r = 5` for :math:`r\in[20,200]` and for the correlation function we'll use a linear binning file with :math:`\Delta r = 1` for :math:`r\in[0,200]`. **NB**: The correlation function binning file must extend down to :math:`r = 0`::

    python python/write_binning_file_linear.py 36 20 200 radial_binning_cov.csv
    python python/write_binning_file_linear.py 200 0 200 radial_binning_corr.csv
    
(See :ref:`write-binning-file` for documentation).

Let's take a look at these::

    20.00000000     25.00000000
    25.00000000     30.00000000
    30.00000000     35.00000000
    35.00000000     40.00000000
    40.00000000     45.00000000
    45.00000000     50.00000000
    50.00000000     55.00000000
    55.00000000     60.00000000
    ....
    
This all looks good.
