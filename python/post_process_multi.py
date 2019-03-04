## Script to post-process the multi-field integrals computed by the C++ code. This computes the shot-noise rescaling parameters, alpha_i, from data derived covariance matrices.
## We output the data and theory jackknife covariance matrices, in addition to full theory covariance matrices and (quadratic-bias corrected) precision matrices. 
## The effective number of samples, N_eff, is also computed.

import numpy as np
import sys,os

# PARAMETERS
if len(sys.argv)!=9:
    print("Usage: python post_process_multi.py {XI_JACKKNIFE_FILE_11} {XI_JACKKNIFE_FILE_12} {XI_JACKKNIFE_FILE_22} {WEIGHTS_DIR} {COVARIANCE_DIR} {N_MU_BINS} {N_SUBSAMPLES} {OUTPUT_DIR}")
    sys.exit()
        
jackknife_file_11 = str(sys.argv[1])
jackknife_file_12 = str(sys.argv[2])
jackknife_file_22 = str(sys.argv[3])
weight_dir = str(sys.argv[4])
file_root = str(sys.argv[5])
m = int(sys.argv[6])
n_samples = int(sys.argv[7])
outdir = str(sys.argv[8])

# Create output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

## First load jackknife xi estimates from data:
print("Loading correlation function jackknife estimates")
xi_jack_11 = np.loadtxt(jackknife_file_11,skiprows=2)
xi_jack_12 = np.loadtxt(jackknife_file_12,skiprows=2)
xi_jack_22 = np.loadtxt(jackknife_file_22,skiprows=2)
assert xi_jack_11.shape==xi_jack_22.shape==xi_jack_12.shape,'Must have the same number of jackknifes for each field.'

n_bins = xi_jack_11.shape[1] # total bins
n_jack = xi_jack_11.shape[0] # total jackknives
n = n_bins//m # radial bins
xi_jack_all = [xi_jack_11,xi_jack_22]

# First exclude any dodgy jackknife regions
good_jk=np.unique(np.where(np.isfinite(xi_jack_11)&np.isfinite(xi_jack_22)&np.isfinite(xi_jack_12))[0])
print("Using %d out of %d jackknives"%(len(good_jk),n_jack))

# Initialize full data covariance matrix
data_cov = np.zeros([2,2,2,2,n_bins,n_bins])

xi_all = np.zeros([2,2,len(good_jk),n_bins])
weights_all = np.zeros_like(xi_all)

# Load in all xi functions
xi_all[0,0]=xi_jack_11[good_jk]
xi_all[0,1]=xi_jack_12[good_jk]
xi_all[1,0]=xi_jack_12[good_jk]
xi_all[1,1]=xi_jack_22[good_jk]

# Load in all weights:
weight_file11 = weight_dir+'jackknife_weights_n%d_m%d_j%d_11.dat'%(n,m,n_jack)
weight_file12 = weight_dir+'jackknife_weights_n%d_m%d_j%d_12.dat'%(n,m,n_jack)
weight_file22 = weight_dir+'jackknife_weights_n%d_m%d_j%d_22.dat'%(n,m,n_jack)
weights11 = np.loadtxt(weight_file11)[:,1:]
weights12 = np.loadtxt(weight_file12)[:,1:]
weights22 = np.loadtxt(weight_file22)[:,1:]
weights_all[0,0]=weights11[good_jk]
weights_all[0,1]=weights12[good_jk]
weights_all[1,0]=weights12[good_jk]
weights_all[1,1]=weights22[good_jk]

# Compute covariance terms:
tmp_cov = np.zeros([2,2,len(good_jk),n_bins])
for j1 in range(2):
    for j2 in range(2):
        mean_xi = np.sum(xi_all[j1,j2]*weights_all[j1,j2],axis=0)/np.sum(weights_all[j1,j2],axis=0)
        tmp_cov[j1,j2] = weights_all[j1,j2]*(xi_all[j1,j2]-mean_xi)

print("Computing all data covariance matrices")
# Now compute covariance matrix:
for j1 in range(2):
    for j2 in range(2):
        for j3 in range(2):
            for j4 in range(2):
                num = np.matmul(tmp_cov[j1,j2].T,tmp_cov[j3,j4])
                denom = np.matmul(weights_all[j1,j2].T,weights_all[j3,j4])
                data_cov[j1,j2,j3,j4] = num/(np.ones_like(denom)-denom)

def load_matrices(index,field,jack=True):
    """Load intermediate or full autocovariance matrices. 
    The field parameter controls which field covariance matrix to load"""
    if jack:
        cov_root = file_root+'CovMatricesJack/'
    else:
        cov_root = file_root+'CovMatricesAll/'
    suffix2 = '_n%d_m%d_%s%s_%s.txt'%(n,m,field,field,index)
    suffix3 = '_n%d_m%d_%s,%s%s_%s.txt'%(n,m,field,field,field,index)
    suffix4 = '_n%d_m%d_%s%s,%s%s_%s.txt'%(n,m,field,field,field,field,index)
    c2 = np.diag(np.loadtxt(cov_root+'c2'+suffix2))
    c3 = np.loadtxt(cov_root+'c3'+suffix3)
    c4 = np.loadtxt(cov_root+'c4'+suffix4)
    if jack:
        EEaA1 = np.loadtxt(cov_root+'EE1'+suffix2)
        EEaA2 = np.loadtxt(cov_root+'EE2'+suffix2)
        RRaA1 = np.loadtxt(cov_root+'RR1'+suffix2)
        RRaA2 = np.loadtxt(cov_root+'RR2'+suffix2)
    
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

alpha_best = np.zeros(2)

indices = ['11','22']

## Optimize for alpha_1 and alpha_2 separately.
for i,index in enumerate(indices):

    RR_file = weight_dir+'binned_pair_counts_n%d_m%d_j%d_%s.dat'%(n,m,n_jack,index)
    #print("Loading weights file from %s"%weight_file)
    #weights = np.loadtxt(weight_file)[:,1:]
    print("DO WE NEED THESE?")
    print("Loading weights file from %s"%RR_file)
    RR=np.loadtxt(RR_file)

    # Filter out bad jackknifes
    xi_jack = xi_jack_all[i][good_jk]
    #weights = weights[good_jk]

    # Read in data covariance matrix
    this_data_cov = data_cov[i,i,i,i]
    
    # Load in full jackknife theoretical matrices
    print("Loading best estimate of jackknife covariance matrix for field %d"%(i+1))
    c2,c3,c4=load_matrices('full',i+1)

    # Check matrix convergence
    from np.linalg import eigvalsh
    eig_c4 = eigvalsh(c4)
    eig_c2 = eigvalsh(c2)
    if min(eig_c4)<-1.*min(eig_c2):
        print("Jackknife 4-point covariance matrix has not converged properly via the eigenvalue test. Exiting")
        sys.exit()
    
    # Load in partial jackknife theoretical matrices
    c2s,c3s,c4s=[],[],[]
    for j in range(n_samples):
        print("Loading field %d jackknife subsample %d of %d"%(i+1,j+1,n_samples))
        c2,c3,c4=load_matrices(j,i+1)
        c2s.append(c2)
        c3s.append(c3)
        c4s.append(c4)

    # Compute inverted matrix
    def Psi(alpha):
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
        return Psi

    def neg_log_L1(alpha):
        """Return negative log L1 likelihood between data and theory covariance matrices"""
        Psi_alpha = Psi(alpha)
        logdet = np.linalg.slogdet(Psi_alpha)
        if logdet[0]<0:
            # Remove any dodgy inversions
            return np.inf        
        return np.trace(np.matmul(Psi_alpha,this_data_cov))-logdet[1]

    # Now optimize for shot-noise rescaling parameter alpha
    print("Optimizing for the shot-noise rescaling parameter alpha_%d"%(i+1))
    from scipy.optimize import fmin
    optimal_alpha = fmin(neg_log_L1,1.)
    print("Optimization complete for field %d - optimal rescaling parameter is alpha_%d = %.6f"%(i+1,i+1,optimal_alpha))
    
    alpha_best[i]=optimal_alpha

# input indices
I1 = [1,1,1,1,2,2]
I2 = [1,2,2,1,1,2]
I3 = [1,1,2,2,2,2]
I4 = [1,1,1,2,2,2]

def matrix_readin(suffix='full'):
    """Read in multi-field covariance matrices. This returns lists of full and jackknife covariance matrices"""

    ## Define arrays for covariance matrices
    c2s,c2js=[np.zeros([2,2,n*m,n*m]) for _ in range(2)]
    RRs=np.zeros([2,2,n*m])
    diff1s,diff2s,JK_weights=[np.zeros([2,2,n_jack,n*m]) for _ in range(3)]
    c3s,c3js=[np.zeros([2,2,2,n*m,n*m]) for _ in range(2)]
    c4s,c4js=[np.zeros([2,2,2,2,n*m,n*m]) for _ in range(2)]  

    for ii in range(len(I1)):
        index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
        index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
        index2="%d%d"%(I1[ii],I2[ii])
        
        j1,j2,j3,j4=I1[ii]-1,I2[ii]-1,I3[ii]-1,I4[ii]-1 # internal indexing

        # Define input files
        file_root_all=file_root+'CovMatricesAll/'
        file_root_jack=file_root+'CovMatricesJack/'
        rr_true_file =weight_dir+'binned_pair_counts_n%d_m%d_j%d_%s.dat'%(n,m,n_jack,index2)
        weights_file = weight_dir+'jackknife_weights_n%d_m%d_j%d_%s.dat'%(n,m,n_jack,index2)
        
        if suffix=='full':
            counts_file = file_root_all+'total_counts_n%d_m%d_%s.txt'%(n,m,index4)
            # Load total number of counts
            total_counts=np.loadtxt(counts_file)
            print("Reading in integral components for C_{%s}, which used %.2e pairs, %.2e triples and %.2e quads of particles"%(index4,total_counts[0],total_counts[1],total_counts[2]))
        else: 
            print("Reading in integral components for C_{%s}, iteration %s"%(index4,suffix))
            
        # Load jackknife weights
        weights=np.loadtxt(weights_file)[:,1:] 

        # Load pair counts
        rr_true = np.loadtxt(rr_true_file)

        # Load full integrals
        c2=np.diag(np.loadtxt(file_root_all+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix)))
        c3=np.loadtxt(file_root_all+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,suffix))
        c4=np.loadtxt(file_root_all+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,suffix))

        # Load jackknife integrals
        c2j=np.diag(np.loadtxt(file_root_jack+'c2_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix)))
        c3j=np.loadtxt(file_root_jack+'c3_n%d_m%d_%s_%s.txt' %(n,m,index3,suffix))
        c4j=np.loadtxt(file_root_jack+'c4_n%d_m%d_%s_%s.txt' %(n,m,index4,suffix))

        # Define cxj components
        EEaA1 = np.loadtxt(file_root_jack+'EE1_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix))
        EEaA2 = np.loadtxt(file_root_jack+'EE2_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix))
        RRaA1 = np.loadtxt(file_root_jack+'RR1_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix))
        RRaA2 = np.loadtxt(file_root_jack+'RR2_n%d_m%d_%s_%s.txt' %(n,m,index2,suffix))
        w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
        w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
        EEa1 = np.sum(EEaA1,axis=0)
        EEa2 = np.sum(EEaA2,axis=0)
        diff1 = EEaA1-w_aA1*EEa1
        diff2 = EEaA2-w_aA2*EEa2

        # Now save components
        RRs[j1,j2]=rr_true
        JK_weights[j1,j2]=weights
        c2s[j1,j2]=c2
        c2js[j1,j2]=c2j
        diff1s[j1,j2]=diff1
        diff2s[j1,j2]=diff2
        c3s[j2,j1,j3]=c3
        c3js[j2,j1,j3]=c3j
        c4s[j1,j2,j3,j4]=c4
        c4js[j1,j2,j3,j4]=c4j
        
        # Add symmetries:
        if j1!=j3:
            c3s[j2,j3,j1]=c3
            c3js[j2,j3,j1]=c3j
        if j1!=j2:
            c4s[j2,j1,j3,j4]=c4
            c4js[j2,j1,j3,j4]=c4j
        if j3!=j4:
            c4s[j1,j2,j4,j3]=c4
            c4js[j1,j2,j4,j3]=c4j
            if j1!=j2:
                c4s[j2,j1,j4,j3]=c4
                c4js[j2,j1,j4,j3]=c4j
        if ((j1!=j3)or(j2!=j4)):
            c4s[j3,j4,j1,j2]=c4.T
            c4js[j3,j4,j1,j2]=c4j.T
            if j3!=j4:
                c4s[j4,j3,j1,j2]=c4.T
                c4js[j4,j3,j1,j2]=c4j.T
            if j1!=j2:
                c4s[j3,j4,j2,j1]=c4.T
                c4js[j3,j4,j2,j1]=c4j.T
                if j3!=j4:
                    c4s[j4,j3,j2,j1]=c4.T
                    c4js[j4,j3,j2,j1]=c4j.T
                    
    def construct_fields(j1,j2,j3,j4,alpha1,alpha2):
        # Reconstruct the full field for given input fields and rescaling parameters

        # Create kronecker deltas
        d_xw=(j1==j4)
        d_xz=(j1==j3)
        d_yw=(j2==j4)
        d_yz=(j2==j3)

        # Compute disconnected piece
        t1=np.matmul(diff1s[j1,j2].T,diff2s[j3,j4])
        t2=np.asarray(np.matmul(np.asmatrix(RRs[j1,j2]).T,np.asmatrix(RRs[j3,j4])))
        t3=1.-np.matmul(JK_weights[j1,j2].T,JK_weights[j3,j4])
        cxj=t1/(t2*t3)

        full=c4s[j1,j2,j3,j4]+0.25*alpha1*(d_xw*c3s[j1,j2,j3]+d_xz*c3s[j1,j2,j4])+0.25*alpha2*(d_yw*c3s[j2,j1,j3]+d_yz*c3s[j2,j1,j4])+0.5*alpha1*alpha2*(d_xw*d_yz+d_xz*d_yw)*c2s[j1,j2]
        jack=c4js[j1,j2,j3,j4]+0.25*alpha1*(d_xw*c3js[j1,j2,j3]+d_xz*c3js[j1,j2,j4])+0.25*alpha2*(d_yw*c3js[j2,j1,j3]+d_yz*c3js[j2,j1,j4])+0.5*alpha1*alpha2*(d_xw*d_yz+d_xz*d_yw)*c2js[j1,j2]+cxj
        return full,jack

    c_tot = np.zeros([2,2,2,2,n*m,n*m])
    cj_tot = np.zeros([2,2,2,2,n*m,n*m])

    for j1 in range(2):
        for j2 in range(2):
            for j3 in range(2):
                for j4 in range(2):
                    c_tot[j1,j2,j3,j4],cj_tot[j1,j2,j3,j4]=construct_fields(j1,j2,j3,j4,alpha_best[0],alpha_best[1])

    return c_tot,cj_tot

# Load full matrices
c_tot,cj_tot = matrix_readin()

# Load subsampled matrices
c_subsamples,cj_subsamples=[],[]
for i in range(n_samples):
    tmp=matrix_readin(i)
    c_subsamples.append(tmp[0])
    cj_subsamples.append(tmp[1])

# Now compute all precision matrices
prec_tot = np.zeros_like(c_tot)
prec_j_tot = np.zeros_like(cj_tot)
iden = np.eye(len(c_tot[0,0,0,0]))

N_eff = np.zeros([2,2,2,2])
D_est = np.zeros_like(c_tot)

def compute_precision(entire_matrix,subsamples):
    summ=0.
    for i in range(n_samples):
        c_excl_i = np.mean(subsamples[:i]+subsamples[i+1:],axis=0)
        summ+=np.matmul(np.linalg.inv(c_excl_i),subsamples[i])
    D_est = (summ/n_samples-iden)*(n_samples-1.)/n_samples
    logdetD = np.linalg.slogdet(D_est)
    D_value = logdetD[0]*np.exp(logdetD[1]/n_bins)
    N_eff_D = (n_bins+1.)/D_value+1.
    precision = np.matmul(iden-D_est,np.linalg.inv(entire_matrix))
    return precision,N_eff_D,D_est

print("Computing precision matrices and N_eff")
for j1 in range(2):
        for j2 in range(2):
            for j3 in range(2):
                for j4 in range(2):
                    full_subsamples = [c_sub[j1,j2,j3,j4] for c_sub in c_subsamples]
                    jack_subsamples = [cj_sub[j1,j2,j3,j4] for cj_sub in cj_subsamples]
                    prec_tot[j1,j2,j3,j4],N_eff[j1,j2,j3,j4],D_est[j1,j2,j3,j4]=compute_precision(c_tot[j1,j2,j3,j4],full_subsamples)
                    prec_j_tot[j1,j2,j3,j4],_,_=compute_precision(cj_tot[j1,j2,j3,j4],jack_subsamples)    

output_name = outdir+'Rescaled_Multi_Field_Covariance_Matrices_n%d_m%d_j%d.npz'%(n,m,n_jack)
np.savez(output_name,jackknife_theory_covariance=cj_tot,full_theory_covariance=c_tot,
         jackknife_data_covariance=data_cov,shot_noise_rescaling=alpha_best,
         jackknife_theory_precision=prec_j_tot,full_theory_precision=prec_tot,
         N_eff=N_eff, full_theory_D_matrix = D_est,
         individual_theory_covariances = c_subsamples)

print("Saved output covariance matrices as %s"%output_name)
