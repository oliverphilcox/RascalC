# Function to reconstruct the covariance matrix integrals computed by the C++ code.

import numpy as np
import sys

if len(sys.argv)<7:
    print "Please specify input parameters in the form {ROOT_DIR} {R_BINS} {MU_BINS} {N_JACK} {MULTI_FIELD} {SHOT_NOISE_PARAM_1} [{SHOT_NOISE_PARAM_2}]."
    sys.exit()
    
root_dir = str(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
J = int(sys.argv[4])
multi_field=int(sys.argv[5])
a = float(sys.argv[6])

if multi_field==0:
    print("Running for a single field")
elif multi_field==1:
    print("Running for two fields")
    if len(sys.argv)<8:
        print("Please specify shot noise rescaling parameter for the second field")
        sys.exit();
    else:
        a2=float(sys.argv[7])
else:
    print("Please specify 0 or 1 for multi_field parameter")
    sys.exit()
    

if multi_field==0:
    I1,I2,I3,I4=[1],[1],[1],[1]
else:
    I1 = [1,1,1,1,2,2]
    I2 = [1,2,2,1,1,2]
    I3 = [1,1,2,2,2,2]
    I4 = [1,1,1,2,2,2]
    
## Define arrays for covariance matrices
c2s,c2js=[np.zeros([2,2,n*m,n*m]) for _ in range(2)]
RRs=np.zeros([2,2,n*m])
diff1s,diff2s,JK_weights=[np.zeros([2,2,J,n*m]) for _ in range(3)]
c3s,c3js=[np.zeros([2,2,2,n*m,n*m]) for _ in range(2)]
c4s,c4js=[np.zeros([2,2,2,2,n*m,n*m]) for _ in range(2)]
  
raise Exception("Remove symmetry factors for single-field case")
  
for ii in range(len(I1)):
    index4="%d%d,%d%d"%(I1[ii],I2[ii],I3[ii],I4[ii])
    index3="%d,%d%d"%(I2[ii],I1[ii],I3[ii])
    index2="%d%d"%(I1[ii],I2[ii])
    j1,j2,j3,j4=I1[ii]-1,I2[ii]-1,I3[ii]-1,I4[ii]-1 # internal indexing
    
    # Define input files
    file_root_all=root_dir+'CovMatricesAll/'
    file_root_jack=root_dir+'CovMatricesJack/'
    rr_true_file =root_dir+'weight_files/binned_pair_counts_n%d_m%d_j%d.dat'%(n,m,J)
    rr_file=file_root_all+'RR_n%d_m%d_%s_full.txt'%(n,m,index2)
    weights_file = root_dir+'weight_files/jackknife_weights_n%d_m%d_j%d.dat'%(n,m,J)
    counts_file = root_dir+'CovMatricesAll/total_counts_n%d_m%d_%s.txt'%(n,m,index4)

    # Load total number of counts
    total_counts=np.loadtxt(counts_file)
    print("Reading in integral components for C_{%s}, which used %.2e pairs, %.2e triples and %.2e quads of particles"%(index4,total_counts[0],total_counts[1],total_counts[2]))
    
    # Load jackknife weights
    weights=np.loadtxt(weights_file)[:,1:] 

    # Load pair counts
    #TODO: Remove rr_est?
    rr_est = np.loadtxt(rr_file)
    rr_true = np.loadtxt(rr_true_file)

    # Load bin counts
    #TODO: Remove?
    binct2=np.loadtxt(file_root_all+'binct_c2_n%d_m%d_%s_full.txt'%(n,m,index2))
    binct3=np.loadtxt(file_root_all+'binct_c3_n%d_m%d_%s_full.txt'%(n,m,index3))
    binct4=np.loadtxt(file_root_all+'binct_c4_n%d_m%d_%s_full.txt'%(n,m,index4))
    binct3=(binct3+binct3.T)/2.
    binct4=(binct4+binct4.T)/2.

    # Load full integrals
    c2=np.diag(np.loadtxt(file_root_jack+'c2_n%d_m%d_%s_full.txt' %(n,m,index2)))
    c3=np.loadtxt(file_root_all+'c3_n%d_m%d_%s_full.txt' %(n,m,index3))
    c4=np.loadtxt(file_root_all+'c4_n%d_m%d_%s_full.txt' %(n,m,index4))
    errc4=np.loadtxt(file_root_all+'c4err_n%d_m%d_%s_full.txt' %(n,m,index4))

    # Load jackknife integrals
    c2j=np.diag(np.loadtxt(file_root_jack+'c2_n%d_m%d_%s_full.txt' %(n,m,index2)))
    c3j=np.loadtxt(file_root_jack+'c3_n%d_m%d_%s_full.txt' %(n,m,index3))
    c4j=np.loadtxt(file_root_jack+'c4_n%d_m%d_%s_full.txt' %(n,m,index4))
    errc4j=np.loadtxt(file_root_jack+'c4err_n%d_m%d_%s_full.txt' %(n,m,index4))
    
    # Define cxj components
    EEaA1 = np.loadtxt(file_root_jack+'EE1_n%d_m%d_%s_full.txt' %(n,m,index2))
    EEaA2 = np.loadtxt(file_root_jack+'EE2_n%d_m%d_%s_full.txt' %(n,m,index2))
    RRaA1 = np.loadtxt(file_root_jack+'RR1_n%d_m%d_%s_full.txt' %(n,m,index2))
    RRaA2 = np.loadtxt(file_root_jack+'RR2_n%d_m%d_%s_full.txt' %(n,m,index2))
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
    
# Reconstruct single field integrals
if multi_field==0:
    c2=c2s[0,0]
    c2j=c2js[0,0]
    diff1=diff1s[0,0]
    diff2=diff2s[0,0]
    c3=c3s[0,0,0]
    c3j=c3js[0,0,0]
    c4=c4s[0,0,0,0]
    c4j=c4js[0,0,0,0]
    RR=RRs[0,0]
    JK_weight=JK_weights[0,0]
    RRaRRb=np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
    fact=np.eye(n*m)-np.matmul(np.asmatrix(JK_weight).T,np.asmatrix(JK_weight)) # 1 - SUM_A(w_aA*w_bA)
    cxj = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
    
    # Define total field
    c_tot = c4+a*c3+a**2.*c2
    cj_tot = c4j+cxj+a*c3j+a**2.*c2j

    # Take transpose symmetry
    #TODO: Remove components?
    c3=(c3+c3.T)/2.
    c2=(c4+c4.T)/2.
    c3j=(c3j+c3j.T)/2.
    c4j=(c4j+c4j.T)/2.
    c_tot=(c_tot+c_tot.T)/2.
    cj_tot=(cj_tot+cj_tot.T)/2.
    
    np.savez(root_dir+"SingleFieldProcessedIntegrals.npz",c_tot=c_tot,cj_tot=cj_tot,shot_noise_rescaling=a)
    print("Output matrices written to %sSingleFieldProcessedIntegrals.npz" %root_dir)
    
    
else:
    print("\nReconstructing multi-field integrals")
    ## Reconstruct general case
    def construct_fields(j1,j2,j3,j4,a,a2):
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
        
        full=c4s[j1,j2,j3,j4]+0.25*a*(d_xw*c3s[j1,j2,j3]+d_xz*c3s[j1,j2,j4])+0.25*a2*(d_yw*c3s[j2,j1,j3]+d_yz*c3s[j2,j1,j4])+0.5*a*a2*(d_xw*d_yz+d_xz*d_yw)*c2s[j1,j2]
        jack=c4js[j1,j2,j3,j4]+0.25*a*(d_xw*c3js[j1,j2,j3]+d_xz*c3js[j1,j2,j4])+0.25*a2*(d_yw*c3js[j2,j1,j3]+d_yz*c3js[j2,j1,j4])+0.5*a*a2*(d_xw*d_yz+d_xz*d_yw)*c2js[j1,j2]+cxj
        return full,jack
    
    c_tot = np.zeros([2,2,2,2,n*m,n*m])
    cj_tot = np.zeros([2,2,2,2,n*m,n*m])
    
    for j1 in range(2):
        for j2 in range(2):
            for j3 in range(2):
                for j4 in range(2):
                    c_tot[j1,j2,j3,j4],cj_tot[j1,j2,j3,j4]=construct_fields(j1,j2,j3,j4,a,a2)
    
        
    np.savez(root_dir+"MultiFieldProcessedIntegrals.npz",c_tot=c_tot,cj_tot=cj_tot,shot_noise_rescaling1=a,shot_noise_rescaling2=a2)
    print("\nOutput matrices writeen to %sMultiFieldProcessedIntegrals.npz"%root_dir)
