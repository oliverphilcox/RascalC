# script to compute covariance matrices for each QPM mock using the relevant correlation function. Also compute for the mean correlation function. This runs the COMAJE C++ code 100 times. 

import numpy as np
import subprocess
import time

init = time.time()

# Compile C++ code
cwd = "/home/oliverphilcox/COMAJE/"
subprocess.run(["bash","clean"],cwd=cwd)
subprocess.run(["make"],cwd=cwd)

all_N_gal = np.zeros(99)

# First compute all N_gal values:
for index in range(99):
    root_dir_m = '/mnt/store1/oliverphilcox/CMU/QPM for Oliver/QPM_Pairs_Mariana/'
    mfile_norms = root_dir_m+'qpm-unrecon-%s-norm.dat'%str(index).zfill(4)
    norms_M = np.loadtxt(mfile_norms,usecols=1)
    factor = norms_M[0]/norms_M[1]
    
    # Define parameters
    all_N_gal[index] =factor*6500000; 
    print(all_N_gal[index])

# Now run code for the mean of the mocks
N_gal = np.mean(all_N_gal)
xi_file = "/mnt/store1/oliverphilcox/CMU/xi_functions/QPM_Mariana_mean.xi"
outfile = "/mnt/store1/oliverphilcox/CMU/QPM_Covariances/Mean/"
    
print("RUNNING FOR MEAN XI \n\n")
subprocess.run(["./cov","-norm", "%.8f"%N_gal,"-cor",xi_file,"-output",outfile],cwd=cwd)

this_time = time.time()-init
print("\n\nMean of mocks finished after %d seconds (%.2f minutes = %.2f hours)\n\n"%(this_time,this_time/60.,this_time/3600.))

for index in range(99):
    # Iterate over QPM index
    root_dir_m = '/mnt/store1/oliverphilcox/CMU/QPM for Oliver/QPM_Pairs_Mariana/'
    mfile_norms = root_dir_m+'qpm-unrecon-%s-norm.dat'%str(index).zfill(4)
    norms_M = np.loadtxt(mfile_norms,usecols=1)
    factor = norms_M[0]/norms_M[1]
    
    # Define parameters
    N_gal =factor*6500000; 
    xi_file = "/mnt/store1/oliverphilcox/CMU/xi_functions/QPM_Mariana_mock_%d.xi"%index
    outfile = "/mnt/store1/oliverphilcox/CMU/QPM_Covariances/Mock_%d/"%index
    
    # Run the code with the correct parameters:
    print("RUNNING FOR INDEX %d OF 99\n\n"%index)
    subprocess.run(["./cov","-norm", "%.8f"%N_gal,"-cor",xi_file,"-output",outfile],cwd=cwd)
    
    this_time = time.time()-init
    print("\n\Mock index %d finished after %d seconds (%.2f minutes = %.2f hours)\n\n"%(index,this_time,this_time/60.,this_time/3600.))
