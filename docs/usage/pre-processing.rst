Pre-Processing
===============

We provide a suite of Python scripts to create input files for the RascalC code. These are found in the ``python/`` directory.

.. _coord-conversion:

Coordinate Conversion
----------------------
This converts a set of input particles (either random particle or galaxy positions) in the form {RA,Dec,redshift,w} to the Cartesian form {x,y,z,w} (for weights w), given some cosmological parameters (using a WCDM coordinate converter by Daniel Eisenstein). The output coordinates are in comoving Mpc/h units, and are saved as an ASCII file for later use.

**Usage**::

    python python/convert_to_xyz.py {INFILE} {OUTFILE} [{OMEGA_M} {OMEGA_K} {W_DARK_ENERGY}]
    
**Parameters**:

- {INFILE}: Input data file containing {RA,Dec,redshift,weight} coordinates for each particle. This may be specified as a ``.fits``, ``.csv``, ``.txt`` or ``.dat`` datafile. For an input FITS file, the columns 'RA', 'DEC', 'Z' and 'WEIGHT_FKP' are required (as in BOSS DR12 data).
- {OUTFILE}: Output ``.txt``, ``.dat`` or ``.csv`` filename.
- *Optional* {OMEGA_M}: Current matter density, :math:`\Omega_m` (default 0.31)
- *Optional* {OMEGA_K}: Current curvature density. :math:`\Omega_k` (default 0)
- *Optional* {W_DARK_ENERGY}: Dark energy equation of state parameter, :math:`w_\Lambda` (default -1)

.. _create-jackknives:

Adding Jackknives
-----------------
This function assigns each particle (either random particles or galxy positions) to a jackknife region, j, by assigning a HEALPix pixel number to each datapoint, with a given value of NSIDE. Data is saved as an ASCII file with {x,y,z,w,j} columns. 

**Usage**::

    python python/create_jackknives.py {INFILE} {OUTFILE} {HEALPIX_NSIDE}
    
**Parameters**:

- {INFILE}: Input data ASCII file of (random/galaxy) Cartesian particle positions with space-separated columns {x,y,z,w}, such as that created by the :ref:`coord-conversion` script. This can be in ``.txt``, ``.dat`` or ``.csv`` format.
- {OUTFILE}: Output ``.txt``, ``.dat`` or ``.csv`` filename.
- {HEALPIX_NSIDE}: HealPix NSIDE parameter which controls the number of pixels used to divide up the sky. For NSIDE=:math:`n`, a total of :math:`12n^2` pixels are used.

.. _particle-subset:

Take Subset of Particles
-------------------------
A utility function to reduce the number of particles in an input ASCII file to a given number. This is primarily used to select a random subsample of the random particle file to speed computation.

**Usage**::

    python pythom/take_subset_of_particles.py {INFILE} {OUTFILE} {N_PARTICLES}
    
**Parameters**:

- {INFILE}: Input data ASCII file with particle positions, in {x,y,z,w}, {x,y,z,w,j} or {RA,Dec,redshift,w} format.
- {OUTFILE}: Outfile ``.txt``, ``.dat`` or ``.csv`` filename.
- {N_PARTICLES}: Desired number of particles in output file. A random sample of length N_PARTICLES is selected from the input file.

.. _write-binning-file:

Create Binning Files
--------------------
A utility function to create radial binning files used by RascalC. We provide three scripts for different binning regimes. This file may be alternatively specified by the user, in the format described in :ref:`file-inputs`.

- *Linear*: Bins are linearly spaced bins in :math:`(r_\mathrm{min},r_\mathrm{max})`.
- *Logarithmic*: Bins are evenly spaced in logarithmic space (base :math:`e`) in :math:`(r_\mathrm{min},r_\mathrm{max})`.
- *Hybrid*: Bins are logarithmically spaced in :math:`(r_\mathrm{min},r_\mathrm{cutoff})`, then linearly spaced in :math:`(r_\mathrm{cutoff},r_\mathrm{max})`.

**Usage**::

    python python/write_binning_file_linear.py {N_BINS} {MIN_R} {MAX_R} {OUTPUT_FILE}
    python python/write_binning_file_logarithmic.py {N_BINS} {MIN_R} {MAX_R} {OUTPUT_FILE}
    python python/wrtie_binning_file_hybrid.py {N_LOG_BINS} {N_LIN_BINS} {MIN_R} {CUTOFF_R} {MAX_R} {OUTPUT_FILE}
    
**Parameters**:

- {N_BINS}: Total number of linear or logarithmic bins.
- {MIN_R}: Minimum bin radius, :math:`r_\mathrm{min}`.
- {MAX_R}: Maximm bin radius, :math:`r_\mathrm{max}`.
- {N_LOG_BINS}: Number of logarithmic bins for hybrid binning.
- {N_LIN_BINS}: Numer of linear bins for hybrid binning.
- {CUTOFF_R}: Radius at which to switch from logarithmic to linear binning, :math:`r_\mathrm{cutoff}` (for hybrid binning).
- {OUTPUT_FILE}: Output ``.txt``, ``.csv`` or ``.dat`` filename.
