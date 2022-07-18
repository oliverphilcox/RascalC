## Script to take a subset of integrals
## More specifically, copy part of partial results to other directory and recompute totals by averaging
## Determines single-field vs multi-field and jackknife automatically

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv)!=6:
    print("Usage: python take_subset_of_integrals.py {COVARIANCE_INPUT_DIR} {N_R_BINS} {N_MU_BINS} {N_SUBSAMPLES_TO_USE} {COVARIANCE_OUTPUT_DIR}")
    sys.exit()

input_root = str(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
n_samples = int(sys.argv[4])
output_root = str(sys.argv[5])

input_root_all = os.path.join(input_root, 'CovMatricesAll/')
input_root_jack = os.path.join(input_root, 'CovMatricesJack/')

output_root_all = os.path.join(output_root, 'CovMatricesAll/')
output_root_jack = os.path.join(output_root, 'CovMatricesJack/')

# Create output directories
os.makedirs(output_root_all, exist_ok=1)
if os.path.exists(input_root_jack): os.makedirs(output_root_jack, exist_ok=1)

# input indices
I1 = [1,1,1,1,1,2,2]
I2 = [1,2,2,2,1,1,2]
I3 = [1,1,2,1,2,2,2]
I4 = [1,1,1,2,2,2,2]

for ii in range(len(I1)): # loop over all field combinations
    index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
    index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
    index2="%d%d"%(I1[ii],I2[ii])

    # full integrals
    c2, c3, c4 = [], [], []
    # read
    for i in range(n_samples):
        try:
            c2.append(np.loadtxt(input_root_all+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        except FileNotFoundError: break # end loop if c2 full not found
        c3.append(np.loadtxt(input_root_all+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i)))
        c4.append(np.loadtxt(input_root_all+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i)))
    if len(c2) == 0: break # end loop if no full integral has been found
    c2, c3, c4 = np.array(c2), np.array(c3), np.array(c4)
    # write
    for i in range(n_samples):
        np.savetxt(output_root_all+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i), c2[i])
        np.savetxt(output_root_all+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i), c3[i])
        np.savetxt(output_root_all+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i), c4[i])
    # average and save
    c2, c3, c4 = np.mean(c2, axis=0), np.mean(c3, axis=0), np.mean(c4, axis=0)
    np.savetxt(output_root_all+'c2_n%d_m%d_%s_full.txt' %(n,m,index2), c2)
    np.savetxt(output_root_all+'c3_n%d_m%d_%s_full.txt' %(n,m,index3), c3)
    np.savetxt(output_root_all+'c4_n%d_m%d_%s_full.txt' %(n,m,index4), c4)
    print("Done with %s full" % index4)

    # jackknife integrals
    c2j, c3j, c4j = [], [], []
    EEaA1, EEaA2, RRaA1, RRaA2 = [], [], [], []
    # read
    for i in range(n_samples):
        try:
            c2j.append(np.loadtxt(input_root_jack+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        except FileNotFoundError: break # end loop if c2 jack not found
        c3j.append(np.loadtxt(input_root_jack+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i)))
        c4j.append(np.loadtxt(input_root_jack+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i)))
        # cxj components
        EEaA1.append(np.loadtxt(input_root_jack+'EE1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        EEaA2.append(np.loadtxt(input_root_jack+'EE2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        RRaA1.append(np.loadtxt(input_root_jack+'RR1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        RRaA2.append(np.loadtxt(input_root_jack+'RR2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
    if len(c2j) == 0: continue # skip rest of the loop if no jack integral has been found
    c2j, c3j, c4j = np.array(c2j), np.array(c3j), np.array(c4j)
    EEaA1, EEaA2, RRaA1, RRaA2 = np.array(EEaA1), np.array(EEaA2), np.array(RRaA1), np.array(RRaA2)
    # write
    for i in range(n_samples):
        np.savetxt(output_root_jack+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i), c2j[i])
        np.savetxt(output_root_jack+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i), c3j[i])
        np.savetxt(output_root_jack+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i), c4j[i])
        # cxj components
        np.savetxt(output_root_jack+'EE1_n%d_m%d_%s_%s.txt' %(n,m,index2,i), EEaA1[i])
        np.savetxt(output_root_jack+'EE2_n%d_m%d_%s_%s.txt' %(n,m,index2,i), EEaA2[i])
        np.savetxt(output_root_jack+'RR1_n%d_m%d_%s_%s.txt' %(n,m,index2,i), RRaA1[i])
        np.savetxt(output_root_jack+'RR2_n%d_m%d_%s_%s.txt' %(n,m,index2,i), RRaA2[i])
    # average and save
    c2j, c3j, c4j = np.mean(c2j, axis=0), np.mean(c3j, axis=0), np.mean(c4j, axis=0)
    EEaA1, EEaA2, RRaA1, RRaA2 = np.mean(EEaA1, axis=0), np.mean(EEaA2, axis=0), np.mean(RRaA1, axis=0), np.mean(RRaA2, axis=0)
    np.savetxt(output_root_jack+'c2_n%d_m%d_%s_full.txt' %(n,m,index2), c2j)
    np.savetxt(output_root_jack+'c3_n%d_m%d_%s_full.txt' %(n,m,index3), c3j)
    np.savetxt(output_root_jack+'c4_n%d_m%d_%s_full.txt' %(n,m,index4), c4j)
    # cxj components
    np.savetxt(output_root_jack+'EE1_n%d_m%d_%s_full.txt' %(n,m,index2), EEaA1)
    np.savetxt(output_root_jack+'EE2_n%d_m%d_%s_full.txt' %(n,m,index2), EEaA2)
    np.savetxt(output_root_jack+'RR1_n%d_m%d_%s_full.txt' %(n,m,index2), RRaA1)
    np.savetxt(output_root_jack+'RR2_n%d_m%d_%s_full.txt' %(n,m,index2), RRaA2)
    print("Done with %s jack" % index4)