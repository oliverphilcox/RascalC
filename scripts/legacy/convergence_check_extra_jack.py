## Script to perform an extra convergence check on jackknife integrals
## More specifically, divide integral subsamples into halves and check similarity of their average results
## Determines single-field vs multi-field automatically
## Now legacy, jackknife integrals convergence checking implemented in convergence_check_extra.py, but old post-processing results do not have the needed information saved in .npz files

import numpy as np
import sys, os


sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "../..")))
from RascalC.convergence_check_extra import convergence_check_extra_splittings


def symmetrized(A: np.ndarray):
    return 0.5 * (A + A.T)


def convergence_check_extra_jack(n: int, m: int, input_root: str, n_samples: int, weight_file: str, RR_file: str, alpha: float = 1, skip_r_bins: int = 0, print_function = print):
    skip_bins = skip_r_bins * m
    input_root_jack = os.path.join(input_root, 'CovMatricesJack/')

    weights = np.loadtxt(weight_file)[:, 1+skip_bins:]
    RR = np.loadtxt(RR_file)

    # input indices
    I1 = [1,1,1,1,1,2,2]
    I2 = [1,2,2,2,1,1,2]
    I3 = [1,1,2,1,2,2,2]
    I4 = [1,1,1,2,2,2,2]

    for ii in range(len(I1)): # loop over all field combinations
        index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
        index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
        index2="%d%d"%(I1[ii],I2[ii])

        # jackknife integrals
        c2j, c3j, c4j = [], [], []
        EEaA1, EEaA2, RRaA1, RRaA2 = [], [], [], []
        # read
        for i in range(n_samples):
            try:
                c2j.append(np.diag(np.loadtxt(input_root_jack+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,i))))
            except (FileNotFoundError, IOError): break # end loop if c2 jack not found
            c3j.append(np.loadtxt(input_root_jack+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,i)))
            c4j.append(np.loadtxt(input_root_jack+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,i)))
            # cxj components
            EEaA1.append(np.loadtxt(input_root_jack+'EE1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
            EEaA2.append(np.loadtxt(input_root_jack+'EE2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
            RRaA1.append(np.loadtxt(input_root_jack+'RR1_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
            RRaA2.append(np.loadtxt(input_root_jack+'RR2_n%d_m%d_%s_%s.txt' %(n,m,index2,i)))
        if len(c2j) == 0: break # terminate the loop if no jack integral has been found
        if len(c2j) < n_samples:
            print_function("Some %s jack samples missing: expected %d, found %d" % (index4, n_samples, len(c2j)))
            continue # skip the rest of the loop like above
        c2j, c3j, c4j = np.array(c2j), np.array(c3j), np.array(c4j)
        EEaA1, EEaA2, RRaA1, RRaA2 = np.array(EEaA1), np.array(EEaA2), np.array(RRaA1), np.array(RRaA2)
        # construct the disconnected term
        RRaRRb = np.matmul(np.asmatrix(RR).T, np.asmatrix(RR))
        fact = np.ones_like(c4j[0]) - np.matmul(np.asmatrix(weights).T, np.asmatrix(weights))
        for i in range(n_samples):
            w_aA1 = RRaA1[i] / np.sum(RRaA1[i], axis=0)
            w_aA2 = RRaA2[i] / np.sum(RRaA2[i], axis=0)
            diff1 = EEaA1[i] - w_aA1 * EEaA1[i].sum(axis=0)
            diff2 = EEaA2[i] - w_aA2 * EEaA2[i].sum(axis=0)
            cx = np.asarray(np.matmul(diff1.T, diff2) / np.matmul(fact, RRaRRb))
            c4j[i] += cx
            # symmetrize the matrices
            # c2j[i] = symmetrized(c2j[i]) # unnecessary because c2 are diagonal
            c3j[i] = symmetrized(c3j[i])
            c4j[i] = symmetrized(c4j[i])
        # compute the covariance
        cj_samples = c2j * alpha**2 + c3j * alpha + c4j
        print_function("Jack covariance")
        return convergence_check_extra_splittings(cj_samples, n_samples, print_function)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) not in (7, 8, 9): # if too few or parity is wrong
        print("Usage: python convergence_check_extra_jack.py {N_R_BINS} {N_MU_BINS} {COVARIANCE_INPUT_DIR} {N_SUBSAMPLES} {JACKKNIFE_WEIGHTS_FILE} {RR_FILE} [{ALPHA} [{SKIP_R_BINS}]]")
        sys.exit(1)

    n = int(sys.argv[1])
    m = int(sys.argv[2])
    input_root = str(sys.argv[3])
    n_samples = int(sys.argv[4])
    weight_file = str(sys.argv[5])
    RR_file = str(sys.argv[6])
    alpha = float(sys.argv[7]) if len(sys.argv) > 7 else 1
    skip_r_bins = int(sys.argv[8]) if len(sys.argv) > 8 else 0

    convergence_check_extra_jack(n, m, input_root, n_samples, weight_file, RR_file, alpha, skip_r_bins)