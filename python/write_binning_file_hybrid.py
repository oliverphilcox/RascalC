## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys
import numpy as np

if len(sys.argv)<7:
    print "Please specify input parameters in the form {N_LOG_BINS} {N_LIN_BINS} {MIN_R} {CUTOFF_R} {MAX_R} {OUTPUT_FILE}."
    sys.exit()
n_log_bins = int(sys.argv[1])
n_lin_bins = int(sys.argv[2])
r_min = float(sys.argv[3])
r_cut = float(sys.argv[4])
r_max = float(sys.argv[5])
out_file = str(sys.argv[6])

print("Using hybrid binning: log binning up to r = %.1f, linear binning up to r = %.1f" %(r_cut,r_max))

# Define radial bins
rbins = np.concatenate((np.logspace(np.log(r_min),np.log(r_cut),n_log_bins+1,base=np.e),np.linspace(r_cut,r_max,n_lin_bins+1)))

# PRINT binning:
with file(out_file,'w+') as writefile:
    for i in range(n_log_bins):
        writefile.write("%.8f\t%.8f\n" %(rbins[i],rbins[i+1]))
    for i in range(n_log_bins+1,n_lin_bins+n_log_bins):
        writefile.write("%.8f\t%.8f\n" %(rbins[i],rbins[i+1]))
    writefile.write("%.8f\t%.8f" %(rbins[n_lin_bins+n_log_bins],rbins[n_lin_bins+n_log_bins+1]))
print "Binning file '%s' written successfully."%out_file
