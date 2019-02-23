# script to count number of galaxies in each QPM mock + produces an npy output file 

import numpy as np
import sys,os

if len(sys.argv)!=3:
    print("USAGE: python count_galaxies.py MIN_MOCK MAX_MOCK")
    sys.exit()
min_mock = int(sys.argv[1])
max_mock = int(sys.argv[2])

assert(min_mock==1)

N_gals = np.zeros(max_mock-min_mock+1,dtype=int)
for mock_no in range(min_mock,max_mock+1):
    print("Counting number of galaxies for mock %d of %d"%(mock_no,max_mock))
    gal_file = '/mnt/store1/oliverphilcox/DR12_QPM/unprocessed/mock_galaxy_DR12_CMASS_N_QPM_%s.rdzw'%str(mock_no).zfill(4)
    with open(gal_file) as infile:
        for l,line in enumerate(infile):
            pass
    N_gals[mock_no-1]=l+1
    
np.save('/mnt/store1/oliverphilcox/DR12_QPM/N_gal.npy',N_gals)
print("Counting complete")
