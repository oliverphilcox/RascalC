"""Convenience script to convert an input (Ra,Dec,w) FITS or txt file to comoving (x,y,z) coordinates saved as a .txt file for use with the main C++ code. (Oliver Philcox 2018 with modifications by Michael Rashkovetskyi 2022, using Daniel Eisenstein's 2015 WCDM Coordinate Converter).
Output file format has (x,y,z,w) coordinates in Mpc/h units 

    Parameters:
        INFILE = input ASCII or FITS file
        OUTFILE = output .txt or .csv file specifier
        ---OPTIONAL---
        OMEGA_M = Matter density (default 0.31)
        OMEGA_K = Curvature density (default 0.0)
        W_DARK_ENERGY = Dark Energy equation of state parameter (default -1.)
        ---FURTHER OPTIONAL---
        USE_FKP_WEIGHTS = whether to use FKP weights column (default False/0; only applies to (DESI) FITS files)
        MASK = sets bins that all must be set in STATUS for the particle to be selected (default 0, only applies to (DESI) FITS files)
        USE_WEIGHTS = whether to use WEIGHTS column, if not, set unit weights (default True/1)

"""

import numpy as np
from astropy.constants import c as c_light
import astropy.units as u
from .utils import read_particles_fits_file
from .wcdm.wcdm import coorddist


H0_h = 100 * u.km / u.s / u.Mpc # h = 1
D_H_Mpch = (c_light / H0_h).to(u.Mpc).value # Hubble distance in Mpc/h


def convert_to_xyz(ra_dec_z_pos: np.ndarray[float], Omega_m: float = 0.31, Omega_k: float = 0, w_dark_energy: float = -1) -> np.ndarray[float]:
    # default cosmology from the BOSS DR12 2016 clustering paper assuming LCDM
    if ra_dec_z_pos.shape[0] != 3: ra_dec_z_pos = ra_dec_z_pos.T
    if ra_dec_z_pos.shape[0] != 3: raise ValueError("The input position do not appear 3D")

    all_ra, all_dec, all_z = ra_dec_z_pos

    # Compute comoving radius in Hubble distances
    all_comoving_radius = coorddist(all_z, Omega_m, w_dark_energy, Omega_k)

    # Convert to Mpc/h
    comoving_radius_Mpch = all_comoving_radius * D_H_Mpch

    # Convert to polar coordinates in radians
    all_phi_rad = np.deg2rad(all_ra)
    all_theta_rad = np.deg2rad(90 - all_dec)

    # Now convert to x,y,z coordinates
    all_z = comoving_radius_Mpch * np.cos(all_theta_rad)
    all_x = comoving_radius_Mpch * np.sin(all_theta_rad) * np.cos(all_phi_rad)
    all_y = comoving_radius_Mpch * np.sin(all_theta_rad) * np.sin(all_phi_rad)
    return np.array((all_x, all_y, all_z)).T

def convert_to_xyz_files(input_file: str, output_file: str, Omega_m: float = 0.31, Omega_k: float = 0, w_dark_energy: float = -1, FKP_weights: bool | tuple[float, str] = False, mask: int = 0, use_weights: bool = True, print_function = print):
    # default cosmology from the BOSS DR12 2016 clustering paper assuming LCDM
    print_function(f"Using cosmological parameters as Omega_m = {Omega_m}, Omega_k = {Omega_k}, w = {w_dark_energy}")
    # Load in data:
    print_function("Reading input file %s in Ra,Dec,z coordinates\n"%input_file)
    if input_file.endswith(".fits"):
        particles = read_particles_fits_file(input_file, FKP_weights, mask, use_weights)
    else:
        particles = np.loadtxt(input_file)
    # in either case, first index (rows) are particle numbers, second index (columns) are different properties

    print_function("Converting redshift to comoving distances:")
    particles[:, :3] = convert_to_xyz(particles[:, :3], Omega_m, Omega_k, w_dark_energy) # assume first three columns are positions, and the rest should not be touched

    print_function("Writing to file %s:"%output_file)
    np.savetxt(output_file, particles)
    print_function("Output positions (of length %d) written succesfully!" % len(particles))
