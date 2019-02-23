import numpy as np


if __name__=='__main__':
    import tqdm
    import multiprocessing as mp
    
    
    file_root = '/mnt/store1/oliverphilcox/DR12_QPM/output_xifinecorr/mock_5/'
    jackknife_file = '/mnt/store1/oliverphilcox/DR12_QPM/xi_jack/mock_5/xi_jack_corrected_n35_m10_j169_11.dat'
    weight_dir = '/mnt/store1/oliverphilcox/DR12_QPM/'
    m=10
    n_samples = 20

    # script to post-process covariance matrices for each QPM mock + produces an npz output file 
    import sys,os

    xi_jack = np.loadtxt(jackknife_file,skiprows=2)
    n_bins = xi_jack.shape[1] # total bins
    n_jack = xi_jack.shape[0] # total jackknives
    n = n_bins//m # radial bins

    weight_file = weight_dir+'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack)
    RR_file = weight_dir+'binned_pair_counts_n%d_m%d_j%d_11.dat'%(n,m,n_jack)

    print("Loading weights file from %s"%weight_file)
    weights = np.loadtxt(weight_file)[:,1:]

    # First exclude any dodgy jackknife regions
    good_jk=np.unique(np.where(np.isfinite(xi_jack))[0])
    print("Using %d out of %d jackknives"%(len(good_jk),n_jack))
    xi_jack = xi_jack[good_jk]
    weights = weights[good_jk]

    # Compute data covariance matrix
    print("Computing data covariance matrix")
    mean_xi = np.sum(xi_jack*weights,axis=0)/np.sum(weights,axis=0)
    tmp = weights*(xi_jack-mean_xi)
    data_cov = np.matmul(tmp.T,tmp)
    denom = np.matmul(weights.T,weights)
    data_cov /= (np.ones_like(denom)-denom)

    print("Loading weights file from %s"%RR_file)
    RR=np.loadtxt(RR_file)

    def load_matrices(index,jack=True):
        """Load intermediate or full covariance matrices"""
        if jack:
            cov_root = file_root+'CovMatricesJack/'
        else:
            cov_root = file_root+'CovMatricesAll/'
        c2 = np.diag(np.loadtxt(cov_root+'c2_n%d_m%d_11_%s.txt'%(n,m,index)))
        c3 = np.loadtxt(cov_root+'c3_n%d_m%d_1,11_%s.txt'%(n,m,index))
        c4 = np.loadtxt(cov_root+'c4_n%d_m%d_11,11_%s.txt'%(n,m,index))
        if jack:
            EEaA1 = np.loadtxt(cov_root+'EE1_n%d_m%d_11_%s.txt' %(n,m,index))
            EEaA2 = np.loadtxt(cov_root+'EE2_n%d_m%d_11_%s.txt' %(n,m,index))
            RRaA1 = np.loadtxt(cov_root+'RR1_n%d_m%d_11_%s.txt' %(n,m,index))
            RRaA2 = np.loadtxt(cov_root+'RR2_n%d_m%d_11_%s.txt' %(n,m,index))

            # Compute disconnected term
            w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
            w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
            diff1 = EEaA1-w_aA1*EEaA1.sum(axis=0)
            diff2 = EEaA2-w_aA2*EEaA2.sum(axis=0)
            RRaRRb = np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
            fact = np.ones_like(c4)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
            cx = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
            c4+=cx

        # Now symmetrize and return matrices
        return c2,0.5*(c3+c3.T),0.5*(c4+c4.T)

    # Load in full jackknife theoretical matrices
    print("Loading best estimate of jackknife covariance matrix")
    c2,c3,c4=load_matrices('full')

    # Load in partial jackknife theoretical matrices
    c2s,c3s,c4s=[],[],[]
    for i in range(n_samples):
        #print("Loading jackknife subsample %d of %d"%(i+1,n_samples))
        cc2,cc3,cc4=load_matrices(i)
        c2s.append(cc2)
        c3s.append(cc3)
        c4s.append(cc4)

    # Compute inverted matrix
    def Psi(alpha,return_neff=False):
        """Compute precision matrix from covariance matrix, removing quadratic order bias terms."""
        c_tot = c2*alpha**2.+c3*alpha+c4
        partial_cov=[]
        for i in range(n_samples):
            partial_cov.append(alpha**2.*c2s[i]+alpha*c3s[i]+c4s[i])
        tmp=0.
        for i in range(n_samples):
            c_excl_i = np.mean(partial_cov[:i]+partial_cov[i+1:],axis=0)
            tmp+=np.matmul(np.linalg.inv(c_excl_i),partial_cov[i])
        D_est=(n_samples-1.)/n_samples * (-1.*np.eye(n_bins) + tmp/n_samples)
        Psi = np.matmul(np.eye(n_bins)-D_est,np.linalg.inv(c_tot))
        if return_neff:
            slogD = np.linalg.slogdet(D_est)
            D_val = slogD[0]*np.exp(slogD[1]/n_bins)
            n_eff = (n_bins+1.)/D_val+1
            return Psi,n_eff
        return Psi

    def jk_calc(excl):
        print("Using jackknife %d of %d"%(excl+1,len(weights)))
        this_weights = np.asarray(list(weights[:excl])+list(weights[excl+1:]))
        this_weights/=np.sum(this_weights,axis=0)
        this_xi_jack = np.asarray(list(xi_jack[:excl])+list(xi_jack[excl+1:]))
        this_mean_xi = np.sum(this_xi_jack*this_weights,axis=0)/np.sum(this_weights,axis=0)
        this_tmp = this_weights*(this_xi_jack-this_mean_xi)
        this_data_cov = np.matmul(this_tmp.T,this_tmp)
        this_denom = np.matmul(this_weights.T,this_weights)
        this_data_cov /= np.ones_like(this_denom)-this_denom
        
        def this_neg_log_L1(alpha):
            """Return negative log L1 likelihood between data and theory covariance matrices"""
            Psi_alpha = np.linalg.inv(c2*alpha**2.+c3*alpha+c4)#Psi(alpha)
            logdet = np.linalg.slogdet(Psi_alpha)
            if logdet[0]<0:
                # Remove any dodgy inversions
                return np.inf        
            return np.trace(np.matmul(Psi_alpha,this_data_cov))-logdet[1]
        from scipy.optimize import fmin
        print("Optimizing jackknife %d"%(excl+1))
        this_alpha=fmin(this_neg_log_L1,1.)
        w_A = np.mean(weights[excl])
        print("Optimization complete for jackknife %d"%(excl+1))
        return this_alpha,w_A
    
    p=mp.Pool(20)
    output=list(tqdm.tqdm(p.imap(jk_calc,np.arange(169)),total=169))
    
    alpha_jk = [o[0] for o in output]
    w_A_jk = [o[1] for o in output]
    
    np.savez("/mnt/store1/oliverphilcox/DR12_QPM/single_mock_jackknife_alpha_simple.npz",alpha=alpha_jk,weights=w_A_jk)
