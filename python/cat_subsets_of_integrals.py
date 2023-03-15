## Script to take a subset of integrals and/or concatenate samples
## More specifically, copy part of partial results to other directory and recompute totals by averaging
## Determines single-field vs multi-field and jackknife automatically
## Do not use if subsamples have different numbers of pairs/triples/quadruplets

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv)<6: # if too few
    print("Usage: python cat_subsets_of_integrals.py {N_R_BINS} {mN_MU_BINS/lMAX_L} {COVARIANCE_INPUT_DIR1} {N_SUBSAMPLES_TO_USE1} [{COVARIANCE_INPUT_DIR2} {N_SUBSAMPLES_TO_USE2} ...] [{COLLAPSE_FACTOR}] {COVARIANCE_OUTPUT_DIR}")
    sys.exit(1)

n = int(sys.argv[1])
mstr = str(sys.argv[2])
input_roots = [str(s) for s in sys.argv[3:-1:2]]
ns_samples = [int(s) for s in sys.argv[4:-1:2]]
output_root = str(sys.argv[-1])
collapse_factor = int(input_roots.pop()) if len(input_roots) > len(ns_samples) else 1 # recover the collapse factor if present
assert len(ns_samples) == len(output_root), "Number of input dirs and subsamples to use from them must be the same"

n_samples_tot = sum(ns_samples)
assert n_samples_tot % collapse_factor == 0, "Collapse factor must divide the total number of samples"
n_samples_out = n_samples_tot // collapse_factor

input_roots_all = [os.path.join(input_root, 'CovMatricesAll/') for input_root in input_roots]
input_roots_jack = [os.path.join(input_root, 'CovMatricesJack/') for input_root in input_roots]

output_root_all = os.path.join(output_root, 'CovMatricesAll/')
output_root_jack = os.path.join(output_root, 'CovMatricesJack/')

# Create output directories
os.makedirs(output_root_all, exist_ok=1)
if any([os.path.exists(input_root_jack) for input_root_jack in input_roots_jack]): os.makedirs(output_root_jack, exist_ok=1)

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
    for input_root_all, n_samples in zip(input_roots_all, ns_samples):
        for i in range(n_samples):
            try:
                c2.append(np.loadtxt(input_root_all+'c2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
            except (FileNotFoundError, IOError): break # end loop if c2 full not found
            c3.append(np.loadtxt(input_root_all+'c3_n%d_%s_%s_%s.txt' % (n, mstr, index3, i)))
            c4.append(np.loadtxt(input_root_all+'c4_n%d_%s_%s_%s.txt' % (n, mstr, index4, i)))
    if len(c2) == 0: break # end loop if no full integral has been found
    if len(c2) < n_samples_tot:
        print("Some %s full samples missing: expected %d, found %d" % (index4, n_samples_tot, len(c2)))
        break # end loop like above
    c2, c3, c4 = [np.array(a) for a in (c2, c3, c4)]
    if collapse_factor > 1: c2, c3, c4 = [np.mean(a.reshape(-1, collapse_factor, *np.shape(a)[1:]), axis=1) for a in (c2, c3, c4)] # average adjacent chunks of collapse_factor samples
    # write
    for i in range(n_samples_out):
        np.savetxt(output_root_all+'c2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), c2[i])
        np.savetxt(output_root_all+'c3_n%d_%s_%s_%s.txt' % (n, mstr, index3, i), c3[i])
        np.savetxt(output_root_all+'c4_n%d_%s_%s_%s.txt' % (n, mstr, index4, i), c4[i])
    # average and save
    c2, c3, c4 = [np.mean(a, axis=0) for a in (c2, c3, c4)]
    np.savetxt(output_root_all+'c2_n%d_%s_%s_full.txt' % (n, mstr, index2), c2)
    np.savetxt(output_root_all+'c3_n%d_%s_%s_full.txt' % (n, mstr, index3), c3)
    np.savetxt(output_root_all+'c4_n%d_%s_%s_full.txt' % (n, mstr, index4), c4)
    print("Done with %s full" % index4)

    # jackknife integrals
    c2j, c3j, c4j = [], [], []
    EEaA1, EEaA2, RRaA1, RRaA2 = [], [], [], []
    # read
    for input_root_jack, n_samples in zip(input_roots_jack, ns_samples):
        for i in range(n_samples):
            try:
                c2j.append(np.loadtxt(input_root_jack+'c2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
            except (FileNotFoundError, IOError): break # end loop if c2 jack not found
            c3j.append(np.loadtxt(input_root_jack+'c3_n%d_%s_%s_%s.txt' % (n, mstr, index3, i)))
            c4j.append(np.loadtxt(input_root_jack+'c4_n%d_%s_%s_%s.txt' % (n, mstr, index4, i)))
            # cxj components
            EEaA1.append(np.loadtxt(input_root_jack+'EE1_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
            EEaA2.append(np.loadtxt(input_root_jack+'EE2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
            RRaA1.append(np.loadtxt(input_root_jack+'RR1_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
            RRaA2.append(np.loadtxt(input_root_jack+'RR2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i)))
    if len(c2j) == 0: continue # skip rest of the loop if no jack integral has been found
    if len(c2j) < n_samples_tot:
        print("Some %s jack samples missing: expected %d, found %d" % (index4, n_samples_tot, len(c2j)))
        continue # skip the rest of the loop like above
    c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2 = [np.array(a) for a in (c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2)]
    if collapse_factor > 1: c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2 = [np.mean(a.reshape(-1, collapse_factor, *np.shape(a)[1:]), axis=1) for a in (c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2)] # average adjacent chunks of collapse_factor samples
    # write
    for i in range(n_samples_out):
        np.savetxt(output_root_jack+'c2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), c2j[i])
        np.savetxt(output_root_jack+'c3_n%d_%s_%s_%s.txt' % (n, mstr, index3, i), c3j[i])
        np.savetxt(output_root_jack+'c4_n%d_%s_%s_%s.txt' % (n, mstr, index4, i), c4j[i])
        # cxj components
        np.savetxt(output_root_jack+'EE1_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), EEaA1[i])
        np.savetxt(output_root_jack+'EE2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), EEaA2[i])
        np.savetxt(output_root_jack+'RR1_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), RRaA1[i])
        np.savetxt(output_root_jack+'RR2_n%d_%s_%s_%s.txt' % (n, mstr, index2, i), RRaA2[i])
    # average and save
    c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2 = [np.mean(a, axis=0) for a in (c2j, c3j, c4j, EEaA1, EEaA2, RRaA1, RRaA2)]
    np.savetxt(output_root_jack+'c2_n%d_%s_%s_full.txt' % (n, mstr, index2), c2j)
    np.savetxt(output_root_jack+'c3_n%d_%s_%s_full.txt' % (n, mstr, index3), c3j)
    np.savetxt(output_root_jack+'c4_n%d_%s_%s_full.txt' % (n, mstr, index4), c4j)
    # cxj components
    np.savetxt(output_root_jack+'EE1_n%d_%s_%s_full.txt' % (n, mstr, index2), EEaA1)
    np.savetxt(output_root_jack+'EE2_n%d_%s_%s_full.txt' % (n, mstr, index2), EEaA2)
    np.savetxt(output_root_jack+'RR1_n%d_%s_%s_full.txt' % (n, mstr, index2), RRaA1)
    np.savetxt(output_root_jack+'RR2_n%d_%s_%s_full.txt' % (n, mstr, index2), RRaA2)
    print("Done with %s jack" % index4)
