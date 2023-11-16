"""Convenience script to convert an input (Ra,Dec,w) FITS or txt file to comoving (x,y,z) coordinates saved as a .txt file for use with the main C++ code. (Oliver Philcox 2018, using Daniel Eisenstein's 2015 WCDM Coordinate Converter).
Output file format has (x,y,z,w) coordinates in Mpc/h units 

    Parameters:
        INFILE = input ASCII file
        OUTFILE = output .txt or .csv file specifier
        ---OPTIONAL---
        OMEGA_M = Matter density (default 0.31)
        OMEGA_K = Curvature density (default 0.0)
        W_DARK_ENERGY = Dark Energy equation of state parameter (default -1.)

"""

import sys
import numpy as np
from astropy.constants import c as c_light
import astropy.units as u

# Load the wcdm module from Daniel Eisenstein
import os
dirname=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, str(dirname)+'/wcdm/')
import wcdm

def convert_to_xyz(ra_dec_z_pos: np.ndarray[float], Omega_m: float = 0.31, Omega_k: float = 0, w_dark_energy: float = -1) -> np.ndarray[float]:
    # default cosmology from the BOSS DR12 2016 clustering paper assuming LCDM
    if ra_dec_z_pos.shape[0] != 3: ra_dec_z_pos = ra_dec_z_pos.T
    if ra_dec_z_pos.shape[0] != 3: raise ValueError("The input position do not appear 3D")

    all_ra, all_dec, all_z = ra_dec_z_pos

    all_comoving_radius = wcdm.coorddist(all_z, Omega_m, w_dark_energy, Omega_k)

    # Convert to Mpc/h
    H_0_h = 100*u.km/u.s/u.Mpc # to ensure we get output in Mpc/h units
    comoving_radius_Mpc = ((all_comoving_radius/H_0_h*c_light).to(u.Mpc)).value

    # Convert to polar coordinates in radians
    all_phi_rad = np.deg2rad(all_ra)
    all_theta_rad = np.deg2rad(90 - all_dec)

    # Now convert to x,y,z coordinates
    all_z = comoving_radius_Mpc * np.cos(all_theta_rad)
    all_x = comoving_radius_Mpc * np.sin(all_theta_rad) * np.cos(all_phi_rad)
    all_y = comoving_radius_Mpc * np.sin(all_theta_rad) * np.sin(all_phi_rad)
    return np.array((all_x, all_y, all_z)).T

def convert_to_xyz_files(input_file: str, output_file: str, Omega_m: float = 0.31, Omega_k: float = 0, w_dark_energy: float = -1, print_function = print):
    # default cosmology from the BOSS DR12 2016 clustering paper assuming LCDM
    print_function(f"Using cosmological parameters as Omega_m = {Omega_m}, Omega_k = {Omega_k}, w = {w_dark_energy}")
    # Load in data:
    print_function("Reading input file %s in Ra,Dec,z coordinates\n"%input_file)
    particles = np.loadtxt(input_file) # first index (rows) are particle numbers, second index (columns) are different properties

    print_function("Converting redshift to comoving distances:")
    particles[:, :3] = convert_to_xyz(particles[:, :3], Omega_m, Omega_k, w_dark_energy) # assume first three columns are positions, and the rest should not be touched

    print_function("Writing to file %s:"%output_file)
    np.savetxt(output_file, particles)
    print_function("Output positions (of length %d) written succesfully!" % len(particles))

if __name__ == "__main__": # if invoked as a script
    if len(sys.argv) not in (3, 6):
        print("Please specify input arguments in the form convert_to_xyz.py {INFILE} {OUTFILE} [{OMEGA_M} {OMEGA_K} {W_DARK_ENERGY}]")
        sys.exit(1)
            
    # Load file names
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])

    from utils import get_arg_safe
    # Read in optional cosmology parameters
    Omega_m = get_arg_safe(3, float, 0.31)
    Omega_k = get_arg_safe(4, float, 0)
    w_dark_energy = get_arg_safe(5, float, -1)
    # defaults from the BOSS DR12 2016 clustering paper assuming LCDM

    convert_to_xyz_files(input_file, output_file, Omega_m, Omega_k, w_dark_energy)