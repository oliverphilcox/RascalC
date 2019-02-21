# script to compute covariance matrices for each QPM mock using the relevant correlation function. Each iteration runs the RascalC code. 

import numpy as np
import subprocess
import time

init = time.time()

# Compile C++ code
cwd = "/home/oliverphilcox/COMAJE/"

min_QPM = 4
max_QPM = 4

# Load N_gal file
N_gal = np.load('/mnt/store1/oliverphilcox/DR12_QPM/N_gal.npy')
binfile_cf = '/mnt/store1/oliverphilcox/DR12_QPM/radial_binning_corr_fine.csv'
m_cf = 20

# Iterate over QPM index    
for index in range(min_QPM,max_QPM+1):
    # Define parameters
    xi_file_ind = "/mnt/store1/oliverphilcox/DR12_QPM/xi_fine/mock_%d/xi_n90_m20_11.dat"%(index)
    outfile = "/mnt/store1/oliverphilcox/DR12_QPM/output_xifinecorr/mock_%d/"%(index)
    
    # Run the code with the correct parameters:
    print("RUNNING FOR INDEX %d OF %d\n\n"%(index,max_QPM))
    print("Using %d galaxies"%N_gal[index-1])
    subprocess.run(["bash","clean"],cwd=cwd) # clean + recompile for safety
    subprocess.run(["make"],cwd=cwd)
    subprocess.run(["./cov","-cor",xi_file_ind,"-output",outfile,"-norm",str(N_gal[index-1]),"-mbin_cf",str(m_cf),"-binfile_cf",binfile_cf],cwd=cwd)
    
    this_time = time.time()-init
    print("\n\Mock index %d of %d finished after %d seconds (%.2f minutes = %.2f hours)\n\n"%(index,max_QPM,this_time,this_time/60.,this_time/3600.))
