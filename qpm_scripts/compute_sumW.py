import numpy as np
import multiprocessing as mp
import tqdm

def compute_sumW(mock_no):
    print("Computing mock %d"%mock_no)
    mock_file = '/mnt/store1/oliverphilcox/DR12_QPM/unprocessed/mock_galaxy_DR12_CMASS_N_QPM_%s.rdzw'%str(mock_no).zfill(4)
    sum_weight=0.
    with open(mock_file,"r") as infile:
        for l,line in enumerate(infile):
            sum_weight+=float(line.split()[-2])
    return sum_weight

if __name__=='__main__':
    p=mp.Pool(20)
    sumWeight=list(tqdm.tqdm(p.imap(compute_sumW,range(1,1001)),total=1000))
    p.close()
    np.save('/mnt/store1/oliverphilcox/DR12_QPM/sumW_DD.npy',sumWeight)
