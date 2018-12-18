import numpy as np
import matplotlib.pyplot as plt
import cmocean

class CovarianceMatrix:
    """Object holding a covariance matrix for the QPM mocks."""
    def params(self,params):
        attributes=['infile_root','n','m','n_indiv','a','r_bins','weights_file','RR_file']
        for attr in attributes:
            setattr(self,attr,getattr(params,attr))
    def read_all_matrices(self,root='full'):
        """Read in the C_{ab} matrices from file"""
        cov_root = self.infile_root+'CovMatricesAll/'
        c2=cov_root+'c2_n%d_m%d_11_%s.txt'%(self.n,self.m,root)
        c3=cov_root+'c3_n%d_m%d_1,11_%s.txt'%(self.n,self.m,root)
        c4=cov_root+'c4_n%d_m%d_11,11_%s.txt'%(self.n,self.m,root)
        c_tot=np.diag(np.loadtxt(c2))*self.a**2.+np.loadtxt(c3)*self.a+np.loadtxt(c4)
        return (c_tot+c_tot.T)/2.
    def read_jack_matrices(self,root='full'):
        """Read in the jackknife C_^J_{ab} matrices from file"""
        cov_root = self.infile_root+'CovMatricesJack/'
        c2=cov_root+'c2_n%d_m%d_11_%s.txt'%(self.n,self.m,root)
        c3=cov_root+'c3_n%d_m%d_1,11_%s.txt'%(self.n,self.m,root)
        c4=cov_root+'c4_n%d_m%d_11,11_%s.txt'%(self.n,self.m,root)
        c_jack_tot=np.diag(np.loadtxt(c2))*self.a**2.+np.loadtxt(c3)*self.a+np.loadtxt(c4)
        c_jack_tot+=self.compute_disconnected(root=root)
        return (c_jack_tot+c_jack_tot.T)/2.
    def compute_disconnected(self,root='full'):
        """Add in the disconnected xi_{ij}xi_{kl} term"""
        cov_root = self.infile_root+'CovMatricesJack/'
        EEaA1 = np.loadtxt(cov_root+'EE1_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
        EEaA2 = np.loadtxt(cov_root+'EE2_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
        RRaA1 = np.loadtxt(cov_root+'RR1_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
        RRaA2 = np.loadtxt(cov_root+'RR2_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
        w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
        w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
        EEa1 = np.sum(EEaA1,axis=0)
        EEa2 = np.sum(EEaA2,axis=0)
        diff1 = EEaA1-w_aA1*EEa1
        diff2 = EEaA2-w_aA2*EEa2
        weights=np.loadtxt(self.weights_file)[:,1:] 
        RR=np.loadtxt(self.RR_file)
        RRaRRb=np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
        fact=np.eye(self.n*self.m)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
        cxj = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
        return cxj
    def compute_eigenvalues(self):
        """Compute the eigenvalues of the symmetric matrices"""
        self.eigval=np.linalg.eigvalsh(self.c_tot)
        self.eigvalJ = np.linalg.eigvalsh(self.c_jack_tot)
    def reduce_matrix(self):
        """Reduce the matrix using C_{reduced,ab} = C_{ab}/sqrt(C_{aa}C_{bb})"""
        self.reduced_matrix=np.zeros_like(self.c_tot)
        for aa in range(len(self.reduced_matrix)):
            for bb in range(len(self.reduced_matrix[0])):
                self.reduced_matrix[aa,bb]=self.c_tot[aa,bb]/np.sqrt(self.c_tot[aa,aa]*self.c_tot[bb,bb])
    def compute_D_est(self):
        """Compute D_est as defined in O'Connell+18"""
        c_tot_mats=[]
        for i in range(self.n_indiv):
            c_tot_mats.append(self.read_all_matrices(root=str(i)))
        nn = len(c_tot_mats)
        summ=0.
        for i in range(nn):
            c_excl_i = np.mean(c_tot_mats[:i]+c_tot_mats[i+1:],axis=0)
            summ+=np.matmul(np.linalg.inv(c_excl_i),c_tot_mats[i])
        self.D_est = (nn-1.)/nn*(-1.*np.eye(len(self.c_tot))+1./nn*summ)
    def compute_precision(self):
        """Compute the precision matrix, using D_est to get less noisy inversion."""
        if not hasattr(self,'D_est'):
            self.compute_D_est()
        self.prec = np.matmul((np.identity(len(self.D_est))-self.D_est),np.linalg.inv(self.c_tot))
    def compute_N_eff(self):
        """Compute the effective number of mocks from D_est"""
        if not hasattr(self,'D_est'):
            self.compute_D_est()
        slogdetD=np.linalg.slogdet(self.D_est)
        n_bins = len(self.D_est)
        D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
        N_eff_D = (n_bins+1.)/D_value+1.
        print("Total N_eff Estimate: %.6e"%N_eff_D)
        return N_eff_D
    def reduce_precision(self):
        """Reduce the precision matrix using prec_{reduced,ab}=prec_{ab}/sqrt(prec_{aa}prec_{bb})"""
        if not hasattr(self,'prec'):
            self.compute_precision()
        self.red_prec=np.zeros_like(self.c_tot)
        for aa in range(len(self.red_prec)):
            for bb in range(len(self.red_prec[0])):
                self.red_prec[aa,bb]=self.prec[aa,bb]/np.sqrt(self.prec[aa,aa]*self.prec[bb,bb])
    def __init__(self, params):
        self.params(params)
        self.c_tot=self.read_all_matrices()
        self.c_jack_tot=self.read_jack_matrices()
        

class plotting_tools:
    """Class to assist with plotting CovarianceMatrix objects"""
    def __init__(self):
        self.FS=16
    def general_plotter(self,matrix,bin_min=None,bin_max=None,title=None):
        """General plotting tool"""
        vmax=np.percentile(matrix.ravel(),99.9)
        plt.matshow(matrix,vmax=vmax,vmin=-vmax,cmap=cmocean.cm.balance)
        plt.ylabel('Bin ID b',fontsize=self.FS);
        plt.xlabel('Bin ID a',fontsize=self.FS);
        plt.gca().xaxis.tick_bottom()
        plt.xlim([bin_min,bin_max])
        plt.ylim([bin_max,bin_min])
        plt.title(title,fontsize=self.FS+4)
        plt.colorbar();
    def plot_covariance(self,matrix_class,bin_min=None,bin_max=None):
        """Plot a covariance matrix for specified bin ranges"""
        self.general_plotter(matrix_class.c_tot,bin_min=bin_min,bin_max=bin_max,title=r'$C_{ab}$')
    def plot_reduced_covariance(self,matrix_class,bin_min=None,bin_max=None):
        """Plot the reduced covariance matrix"""
        if not hasattr(matrix_class,'reduced_matrix'):
            matrix_class.reduce_matrix()
        self.general_plotter(matrix_class.reduced_matrix,bin_min=bin_min,bin_max=bin_max,title=r'Reduced $C_{ab}$')
    def plot_precision(self,matrix_class,bin_min=None,bin_max=None):
        """Plot and compute the precision matrix"""
        if not hasattr(matrix_class,'prec'):
            matrix_class.compute_precision()
        self.general_plotter(matrix_class.prec,bin_min=bin_min,bin_max=bin_max,title=r'$\psi_{ab}$')
    def plot_reduced_precision(self,matrix_class,bin_min=None,bin_max=None):
        """Plot and compute the reduced precision matrix"""
        if not hasattr(matrix_class,'red_prec'):
            matrix_class.reduce_precision()
        self.general_plotter(matrix_class.red_prec,bin_min,bin_max,title=r'Reduced $\psi_{ab}$')
    def plot_eigenvalues(self,matrix_class):
        """Plot and compute the eigenvalues of the matrices"""
        if not hasattr(matrix_class,'eigval'):
            matrix_class.compute_eigenvalues()
        plt.plot(matrix_class.eigval,label='Full Matrix');
        plt.ylabel('Eigenvalue',fontsize=self.FS);plt.xlabel('Index',fontsize=self.FS)
        plt.plot(matrix_class.eigvalJ,label='Jackknife Matrix')
        plt.legend(fontsize=self.FS-4.);
        plt.yscale('log');
