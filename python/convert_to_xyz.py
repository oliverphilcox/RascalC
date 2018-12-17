"""Convenience script to convert an input (Ra,Dec,w) FITS file to comoving (x,y,z) coordinates saved as a .txt file for use with the main C++ code. (Oliver Philcox 2018, using Daniel Eisenstein's 2015 WCDM Coordinate Converter).
Output file format has (x,y,z,w) coordinates in Mpc/h units 

    Parameters:
        INFILE = input fits/csv/txt/dat file
        OUTFILE = output .txt or .csv file specifier
        FILETYPE = specifies the filetype. Must be 'fits'/'txt'/'csv'/'dat'
        ---OPTIONAL---
        OMEGA_M = Matter density (default 0.31)
        OMEGA_K = Curvature density (default 0.0)
        W_DARK_ENERGY = Dark Energy equation of state parameter (default -1.)

"""

import sys
import numpy as np

# Read in optional cosmology parameters
if len(sys.argv)==7:
    omega_m = float(sys.argv[4])
    omega_k = float(sys.argv[5])
    w_dark_energy = float(sys.argv[6])
elif len(sys.argv)==4: # use defaults (from the BOSS DR12 2016 clustering paper assuming LCDM)
    omega_m = 0.31
    omega_k = 0.
    w_dark_energy = -1.
else:
    print("Please specify input arguments in the form convert_to_xyz.py {INFILE} {OUTFILE} {FILETYPE} [{OMEGA_M} {OMEGA_K} {W_DARK_ENERGY}]")
    sys.exit()

print("\nUsing cosmological parameters as Omega_m = %.2f, Omega_k = %.2f, w = %.2f" %(omega_m,omega_k,w_dark_energy))
          
# Load file names
input_file = str(sys.argv[1])
output_file = str(sys.argv[2])
filetype = str(sys.argv[3])
print("\nUsing input file %s in Ra,Dec,z coordinates\n"%input_file)

# Load the wcdm module from Daniel Eisenstein
import os
dirname=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, str(dirname)+'/wcdm/')
import wcdm

# Load in data:
if filetype=='fits':
    from astropy.io import fits
    fitsfile=fits.open(input_file)[1]
    all_ra=fitsfile.data['RA']
    all_dec=fitsfile.data['DEC']
    all_z=fitsfile.data['Z']
    all_w=fitsfile.data['WEIGHT_FKP']
elif ((filetype=='dat')or(filetype=='txt')or(filetype=='csv')):
    print("Counting lines in file")
    total_lines=0
    for n, line in enumerate(open(input_file, 'r')):
        total_lines+=1

    all_ra,all_dec,all_z,all_w=[np.zeros(total_lines) for _ in range(4)]

    print("Reading in data");
    for n, line in enumerate(open(input_file, 'r')):
        if n%1000000==0:
            print("Reading line %d of %d" %(n,total_lines))
        split_line=np.array(line.strip().split("  "), dtype=float) 
        all_ra[n]=split_line[0];
        all_dec[n]=split_line[1];
        all_z[n]=split_line[2];
        all_w[n]=split_line[3];
else: 
    print("Must specify filetype as 'txt', 'csv', 'dat' or 'fits'")
    sys.exit()
    
from astropy.constants import c as c_light
import astropy.units as u

print("Converting z to comoving distances:")
all_comoving_radius=wcdm.coorddist(all_z,omega_m,w_dark_energy,omega_k)

# Convert to Mpc/h
H_0_h=100*u.km/u.s/u.Mpc # to ensure we get output in Mpc/h units
H_0_SI = H_0_h.to(1./u.s)
comoving_radius_Mpc = ((all_comoving_radius/H_0_SI*c_light).to(u.Mpc)).value

# Convert to polar coordinates in radians
all_phi_rad = all_ra*np.pi/180.
all_theta_rad = np.pi/2.-all_dec*np.pi/180.

# Now convert to x,y,z coordinates
all_z = comoving_radius_Mpc*np.cos(all_theta_rad)
all_x = comoving_radius_Mpc*np.sin(all_theta_rad)*np.cos(all_phi_rad)
all_y = comoving_radius_Mpc*np.sin(all_theta_rad)*np.sin(all_phi_rad)

print("Writing to file %s:"%output_file)
# Now write to file:
with open(output_file,"w+") as outfile:
    for p in range(len(all_z)):
        outfile.write("%.8f %.8f %.8f %.8f\n" %(all_x[p],all_y[p],all_z[p],all_w[p]))
print("Output positions (of length %d) written succesfully!"%len(all_z))
