## Utility function to create the Corrfunc radial-binning file for a given number of radial bins in a certain range of comoving radii.

import sys
import numpy as np

if len(sys.argv)<5:
    print "Please specify input parameters in the form {N_LOG_BINS} {MIN_R} {MAX_R} {OUTPUT_FILE}."
    sys.exit()
nrbins = int(sys.argv[1])
r_min = float(sys.argv[2])
r_max = float(sys.argv[3])
out_file = str(sys.argv[4])

print("Using LOG binning");

# Define radial bins
rbins = np.logspace(np.log(r_min),np.log(r_max),nrbins+1,base=np.e)

# PRINT binning:
with file(out_file,'w+') as writefile:
    for i in range(nrbins):
        writefile.write("%.8f\t%.8f\n" %(rbins[i],rbins[i+1]))
    writefile.write("%.8f\t%.8f" %(rbins[nrbins-1],rbins[nrbins]))
print "Binning file '%s' written successfully."%out_file
