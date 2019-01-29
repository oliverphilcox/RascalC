import numpy as np
import matplotlib.pyplot as plt
import cmocean

class CovarianceMatrix:
    """ This object holds a reconstructed covariance matrix and provides routines for computing useful features such as its eigenvalues and inverses. Both the full and jackknife covariance matrices are provided, and we combine the individual 2-, 3- and 4-point terms.
    
    The covariance matrix parameters can be specified as an instance of the :class:`Parameters` class or as some user defined class with the same attributes.
    
    On class initialization, the attributes :attr:`c_tot` and :attr:`c_jack_tot` are created, which contain the full and jackknife covariance matrices. 
    
    Attributes:
        - :attr:`infile_root` (str): Input file directory (inherited from :class:`Parameters` class)
        - :attr:`n` (int): Number of radial bins (inherited from :class:`Parameters` class)
        - :attr:`m` (int): Number of angular bins (inherited from :class:`Parameters` class)
        - :attr:`n_indiv` (int): Number of individual matrix estimates (inherited from :class:`Parameters` class)
        - :attr:`a` (float): Shot-noise rescaling parameter (inherited from :class:`Parameters` class)
        - :attr:`r_bins` (np.ndarray): Radial binning lower and upper limits (inherited from :class:`Parameters` class)
        - :attr:`weights_file` (str): Jackknife weights filename. Generally of the form ``jackknife_weights_n{}_m{}_j{}.dat`` (inherited from :class:`Parameters` class)
        - :attr:`RR_file` (str): Summed RR pair counts filename. Generally of the form ``RR_pair_counts_n{}_m{}_j{}.dat`` (inherited from :class:`Parameters` class)
        - :attr:`c_tot` (np.ndarray): Full covariance matrix (created on initialization)
        - :attr:`c_jack_tot` (np.ndarray): Jackknife covariance matrix (created on initialization)
        - :attr:`prec` (np.ndarray): Precision matrix (created by :func:`compute_precision`)
        - :attr:`prec_jack` (np.ndarray): Jackknife precision matrix (created by :func:`compute_precision`)
        - :attr:`eigval` (np.ndarray): Real eigenvalues of the total covariance matrix :attr:`c_tot` (created by :func:`compute_eigenvalues`)
        - :attr:`eigvalJ` (np.ndarray): Real eigenvalues of the jackknife covariance matrix :attr:`c_jack_tot` (created by :func:`compute_eigenvalues`)
        - :attr:`reduced_matrix` (np.ndarray): Reduced form of the total covariance matrix :attr:`c_tot` (created by :func:reduce_matrix).
        - :attr:`reduced_jack_matrix` (np.ndarray): Reduced form of the jackknife covariance matrix :attr:`c_jack_tot` (created by :func:reduce_matrix).
        - :attr:`reduced_precision` (np.ndarray): Reduced form of the total covariance matrix :attr:`prec` (created by :func:`reduce_matrix`).
        - :attr:`reduced_precision_jack` (np.ndarray): Reduced form of the jackknife covariance matrix :attr:`prec_jack` (created by :func:`reduce_matrix`)
        - :attr:`D_est` (np.ndarray): :math:`\hat{D}` matrix (computed by :func:`compute_D_est`)
        - :attr:`D_est_jack` (np.ndarray): :math:`\hat{D}` matrix (computed by :func:`compute_D_est_jack`)
    
    Usage:
        
        >>> par = Parameters(infile_root,n,m,n_indiv,a,r_bins,weights_file,RR_file)
        >>> my_matrix = CovarianceMatrix(par)
        >>> my_jackknife_matrix = my_matrix.c_jack_tot
        
    
    """
    def __init__(self, params):
        """ Initialize the :class:`CovarianceMatrix` object and read in full and jackknife matrices from the file locations specified in the :class:`Parameters` object.

        Args:
            - params (:class:`Parameters`):  The :class:`Parameters` object which contains the input file locations and information.

        Returns:
            None
            
        Raises:
            AttributeError
            
        This creates the :attr:`c_tot` and :attr:`c_jack_tot` attributes, which contain the full and jackknife covariance matrices.
    
        """
        
        # First assign parameters
        attributes=['infile_root','n','m','n_indiv','a','r_bins','weights_file','RR_file']
        for attr in attributes:
            try:
                setattr(self,attr,getattr(params,attr))
            except AttributeError:
                raise AttributeError("Attribute %s not found in parameter file"%attr)
        
        # Now read in matrice
        self.c_tot=self.read_matrices(matrix_type='full')
        self.c_jack_tot=self.read_matrices(matrix_type='jack')
        
    def read_matrices(self,root='full',matrix_type='full'):
        """
        Read in :math:`C_{ab}` matrices from file and combine 2-, 3- and 4-point estimates. This reads in the complete file by default. 
        
        .. note:: The combined matrix is read in by default on class initialization
        
        Kwargs:
            - root (str): The matrix index to read in. This can be an integer up to :attr:`CovarianceMatrix.n_indiv`. If root='full', the total combined matrix is read in.
            - matrix_type (str): Either 'full' or 'jack', this specifies what type of matrix to read in (either the total or jackknife matrix).
        
        Returns:
            - cov (np.ndarray): Symmetrized matrix array of dimension (n*m) x (n*m)
        
        Raises:
            FileNotFoundError, Exception
        
        """
        if matrix_type=='full':
            cov_root = self.infile_root+'CovMatricesAll/'
        elif matrix_type=='jack':
            cov_root = self.infile_root+'CovMatricesJack/'
        else:
            raise Exception("matrix_type parameter must be either 'full' or 'jack'")
        
        # Create filenames
        c2=cov_root+'c2_n%d_m%d_11_%s.txt'%(self.n,self.m,root)
        c3=cov_root+'c3_n%d_m%d_1,11_%s.txt'%(self.n,self.m,root)
        c4=cov_root+'c4_n%d_m%d_11,11_%s.txt'%(self.n,self.m,root)
        
        # Read in files
        try:
            cov=np.diag(np.loadtxt(c2))*self.a**2.+np.loadtxt(c3)*self.a+np.loadtxt(c4)
        except FileNotFoundError:
            raise FileNotFoundError("Files with index %s not found in the %s directory"%(cov_root,root))
        
        # Add disconnected term if necessary
        if matrix_type=='jack':
            cov+=self._compute_disconnected(root=root)
        
        # Return symmetrized object
        return (cov+cov.T)/2.
    
    def _compute_disconnected(self,root='full'):
        """ Internal function to compute the disconnected contribution to the jackknife integral.
        This is found from the two independent estimates of the EE and RR parameters found by the C++ code.
        
        Args:
            - root (str): The matrix index to read in. This can be an integer up to :attr:`CovarianceMatrix.n_indiv`. If root='full', the total combined matrix is read in.
            
        Returns:
            - cxj (np.ndarry): Disconnect jackknife term array. 
            
        Raises:
            FileNotFoundError
        """
        cov_root = self.infile_root+'CovMatricesJack/'
        # Load input files
        try:
            EEaA1 = np.loadtxt(cov_root+'EE1_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
            EEaA2 = np.loadtxt(cov_root+'EE2_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
            RRaA1 = np.loadtxt(cov_root+'RR1_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
            RRaA2 = np.loadtxt(cov_root+'RR2_n%d_m%d_11_%s.txt' %(self.n,self.m,root))
        except FileNotFoundError:
            raise FileNotFoundError("Input EEaA and RRaA files not found")
        
        # Compute weight set
        w_aA1 = RRaA1/np.sum(RRaA1,axis=0)
        w_aA2 = RRaA2/np.sum(RRaA2,axis=0)
        EEa1 = np.sum(EEaA1,axis=0)
        EEa2 = np.sum(EEaA2,axis=0)
        diff1 = EEaA1-w_aA1*EEa1
        diff2 = EEaA2-w_aA2*EEa2
        
        # Load jackknife weights
        try:
            weights=np.loadtxt(self.weights_file)[:,1:] 
        except FileNotFoundError:
            raise FileNotFoundError("Weights file %s not found."%self.weights_file)
        try:
            RR=np.loadtxt(self.RR_file)
        except FileNotFoundError:
            raise FileNotFoundError("RR pair count file %s not found."%self.RR)
            
        # Compute disconnected term
        RRaRRb=np.matmul(np.asmatrix(RR).T,np.asmatrix(RR))
        fact=np.eye(self.n*self.m)-np.matmul(np.asmatrix(weights).T,np.asmatrix(weights))
        cxj = np.asarray(np.matmul(diff1.T,diff2)/np.matmul(fact,RRaRRb))
        return cxj
    
    def compute_eigenvalues(self):
        """
        Compute the eigenvalues of the symmetric full and jackknife matrices. These are set to the :attr:`eigval` and :attr:`eigvalJ` attributes and returned. The input matrices are symmetrized so all eigenvalues will be real.
        
        Returns:
            - eigval (np.npdarray): Ordered list of real eigenvalues of the full covariance matrix
            - eigvalJ (np.npdarray): Ordered list of real eigenvalues of the jackknife covariance matrix.
        
        """
        self.eigval=np.linalg.eigvalsh(self.c_tot)
        self.eigvalJ = np.linalg.eigvalsh(self.c_jack_tot)
        
        return self.eigval, self.eigvalJ
    
    def _which_matrix(self,matrix_type):
        """ Utility function to return the relevant matrix given a input type string.
        
        Args:
            - matrix_type (str): Either 'full', 'jack' or 'prec', this specifies what type of matrix to use (either the total, jackknife or precision matrix).
        
        Returns:
            - cov_matrix (np.npdarray): Either full, jackknife or precision covariance matrix.
            
        Raises:
            Exception
        """
        if matrix_type=='full':
            cov_matrix = self.c_tot
        elif matrix_type=='jack':
            cov_matrix = self.c_jack_tot
        elif matrix_type=='prec':
            if not hasattr(self,'prec'):
                self.compute_precision('full')
            cov_matrix = self.prec
        elif matrix_type=='prec_jack':
            if not hasattr(self,'prec_jack'):
                self.compute_precision('jack')
            cov_matrix = self.prec_jack
        else:
            raise Exception("matrix_type parameter must be either 'full', 'jack', 'prec' or 'prec_jack'")
        
        return cov_matrix

    def reduce_matrix(self,matrix_type='full'):
        """
        Function to reduce an input correlation matrix to the form :math:`C_{\mathrm{reduced},ab} = C_{ab}/\sqrt(C_{aa}C_{bb})`. This is returned and set as the relevant attribute; :attr:`reduced_matrix`,  :attr:`reduced_jack_matrix` or :attr:`reduced_precision`.
        
        Args:
            - matrix_type (str): Must be 'full', 'jack', 'prec' or 'prec_jack' specifying whether to reduce the full, jackknife, full precision or jackknife precision matrix.
            
        Returns: 
            - reduced_matrix (np.ndarray): Reduced covariance matrix of the requested type.
            
        Raises:
            Exception
        """
        this_mat = self._which_matrix(matrix_type)
        diag = np.diag(this_mat)
        divisor = np.matmul(diag.reshape((len(diag),1)),diag.reshape((1,len(diag))))
        reduced_matrix = np.divide(this_mat,np.sqrt(divisor))
        if matrix_type=='full':
            self.reduced_matrix = reduced_matrix
        elif matrix_type=='jack':
            self.reduced_jack_matrix = reduced_matrix
        elif matrix_type=='prec':
            self.reduced_precision = reduced_matrix
        elif matrix_type=='prec_jack':
            self.reduced_precision_jack = reduced_matrix
        return reduced_matrix        
        
    def compute_D_est(self,matrix_type='full'):
        """ 
        
        Compute the estimated :math:`\hat{D}` matrix as defined in O'Connell & Eisenstein 2018.
        This is computed as :math:`\frac{ N-1 }{ N }\left(-I +\frac{ 1 }{ N }\sum_i C^{-1}_{[i]}C_i\right)` for the :math:`i`th of :math:`N` covariance matrix estimate :math:`C_i` , where :math:`C_{[i]}` excludes :math:`C_i`. This sets either the :attr:`D_est` or the :attr:`D_est_jack` attribute.
        
        Args:
            - matrix_type (str): Either 'full' or 'jack', this specifies what type of matrix to use (either the total or jackknife matrix).
        
        Returns:
            - D_est (np.ndarray): D matrix estimate. Also set as an attribute.
        
        Raises:
            Exception, FileNotFoundError
        
        
        """
        # Read in individual matrices
        c_mats=[]
        for i in range(self.n_indiv):
            c_mats.append(self.read_matrices(root=str(i),matrix_type=matrix_type))
        nn = len(c_mats)
        
        # Compute summand term
        summ=0.
        for i in range(nn):
            c_excl_i = np.mean(c_mats[:i]+c_mats[i+1:],axis=0)
            summ+=np.matmul(np.linalg.inv(c_excl_i),c_mats[i])
        D_est = (nn-1.)/nn*(-1.*np.eye(len(c_mats[0]))+1./nn*summ)
        if matrix_type=='full':
            self.D_est=D_est
        else:
            self.D_est_jack = D_est
        return D_est

    def compute_precision(self,matrix_type='full'):
        """ Compute the precision matrix, via the method of O'Connell & Eisenstein 2018, using :attr:`D_est` to get less noisy inversion i.e. :math:`\Psi = (I-\hat{D})\times C^{-1}` 
        
        This will compute the :attr:`D_est` matrix if not already present and sets either the :attr:`prec` or :attr:`prec_jack` attribute.
        
        Args:
            - matrix_type (str): Either 'full' or 'jack', this specifies what type of matrix to use (either the total or jackknife matrix).
        
        Returns:
            - prec (np.ndarray): Precision matrix array. 
        
        Raises:
            Exception
        
        """
        if matrix_type=='full':
            if not hasattr(self,'D_est'):
                self.compute_D_est('full')
            prec = np.matmul((np.identity(len(self.D_est))-self.D_est),np.linalg.inv(self.c_tot))
            self.prec = prec
        elif matrix_type=='jack':
            if not hasattr(self,'D_est_jack'):
                self.compute_D_est('jack')
            prec = np.matmul((np.identity(len(self.D_est_jack))-self.D_est_jack),np.linalg.inv(self.c_jack_tot))
            self.prec_jack = prec
        else:
            raise Exception("Matrix type must be 'full' or 'jack'")
        return prec
    
    def compute_N_eff(self,matrix_type='full'):
        """ Compute the effective number of mocks from :attr:`D_est`, using :math:`N_\mathrm{eff} = \frac{n_\mathrm{bins}+1}{|D|^{1/n_\mathrm{bins}} }+1`
        
        Args:
            - matrix_type (str): Either 'full' or 'jack', this specifies what type of matrix to use (either the total or jackknife matrix).
            
        Returns:
            - N_eff_D (float64): Effective number of mocks 
            
        Raises:
            Exception
            
        """
        if matrix_type=='full':
            if not hasattr(self,'D_est'):
                self.compute_D_est('full')
            D_est = self.D_est
        elif matrix_type=='jack':
            if not hasattr(self,'D_est_jack'):
                self.compute_D_est('jack')
            D_est = self.D_est_jack
        else:
            raise Exception("matrix_type must be 'full' or 'jack'")
        slogdetD=np.linalg.slogdet(D_est)
        n_bins = len(D_est)
        D_value = slogdetD[0]*np.exp(slogdetD[1]/n_bins)
        N_eff_D = (n_bins+1.)/D_value+1.
        print("Total N_eff Estimate: %.6e"%N_eff_D)
        return N_eff_D
        
class Parameters:
    """Class containing covariance matrix parameters. The parameters are set on initialization and are used to construct the :class:`CovarianceMatrix` class. 
    
    Attributes:
        - :attr:`infile_root` (str): Input file directory
        - :attr:`n` (int): Number of radial bins
        - :attr:`m` (int): Number of angular bins
        - :attr:`n_indiv` (int): Number of individual matrix estimates
        - :attr:`a` (float): Shot-noise rescaling parameter
        - :attr:`r_bins` (np.ndarray): Radial binning upper and lower bins [array dimension (2xn)]
        - :attr:`weights_file` (str): Jackknife weights filename. Generally of the form ``jackknife_weights_n{}_m{}_j{}.dat``
        - :attr:`RR_file` (str): Summed RR pair counts filename. Generally of the form ``RR_pair_counts_n{}_m{}_j{}.dat``
        
    Usage: 
        >>> par = Parameters(infile_root,n,m,n_indiv,a,r_bins.weight_file,RR_file)
        >>> my_matrix_class = CovarianceMatrix(par)
            
    If multiple covariance matrices are to be read in with similar parameters, we can simply alter attributes in this file e.g.
    
        >>> par.a = 5 # alter this parameter
        >>> my_matrix_class2 = CovarianceMatrix(par) # create new matrix class
    
    """
    def __init__(self,infile_root,n,m,n_indiv,a,r_bins,weights_file,RR_file):
        self.infile_root = infile_root
        self.n = n
        self.m = m
        self.n_indiv = n_indiv
        self.a = a
        self.r_bins = r_bins
        self.weights_file = weights_file
        self.RR_file = RR_file        

class plotting_tools:
    """ Class to assist with plotting :class:`CovarianceMatrix` objects by providing various plotting routines. This uses a :mod:`matplotlib` backend.
    
    Attributes:
        - :attr:`FS` (int): :mod:`matplotlib` fontsize
        
    Usage:
        This class contains various routines for plotting matrix attributes. All are used in the following manner::
            >>> pt = plotting_tools() # load class
            >>> my_matrix = CovarianceMatrix(pars) # load matrix given some parameter class
            >>> plot_covariance(my_matrix) # plot covariance matrix
    """
    def __init__(self,FS=16):
        """ Initialize the class.
        
        Args:
            - FS (int): :mod:`matplotlib` fontsize
        """
        self.FS=FS
    def _general_plotter(self,matrix,bin_min=None,bin_max=None,title=""):
        """General plotting tool for matrices used by other routines in the class. This is based on the :mod:`matplotlib` :func:`matshow` routine.   
        
        Args:
            - matrix (np.ndarray): Matrix to plot
        
        Kwargs:
            - bin_min (int): Minimum bin number
            - bin_max (int): Maximum bin number
            - title (str): Figure title
            
        Returns:
            - figure: :mod:`matplotlib` figure object
        """
        vmax=np.percentile(matrix.ravel(),99.9)
        figure = plt.figure()
        plt.matshow(matrix,vmax=vmax,vmin=-vmax,cmap=cmocean.cm.balance)
        plt.ylabel('Bin ID b',fontsize=self.FS);
        plt.xlabel('Bin ID a',fontsize=self.FS);
        plt.gca().xaxis.tick_bottom()
        plt.xlim([bin_min,bin_max])
        plt.ylim([bin_max,bin_min])
        if len(title)>0:
            plt.title(title,fontsize=self.FS+4)
        plt.colorbar();
        return figure
    def plot_covariance(self,matrix_class,bin_min=None,bin_max=None,matrix_type='full',title=""):
        """ 
        Plot the covariance matrix for specified bin ranges.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object
        
        Kwargs:
            - bin_min (int): Minimum bin number
            - bin_max (int): Maximum bin number
            - matrix_type (str): Matrix type ('full' or 'jack')
            - title (str): Figure title (overriding the default)
        Returns:
            - fig: Matplotlib figure object.
        
        Raises:
            Exception
        """
        if matrix_type=='full':
            mat=matrix_class.c_tot
            name = r'$C_{ab}$'
        elif matrix_type=='jack':
            mat=matrix_class.c_jack_tot
            name = r'$C_{ab}^J$'
        else:
            raise Exception("matrix_type must be 'full' or 'jack'")
        if len(title)>0:
            name = title
        return self._general_plotter(mat,bin_min=bin_min,bin_max=bin_max,title=name)
    def plot_reduced_covariance(self,matrix_class,bin_min=None,bin_max=None):
        """ 
        Plot the reduced covariance matrix for specified bin ranges.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object
        
        Kwargs:
            - bin_min (int): Minimum bin number
            - bin_max (int): Maximum bin number
            - matrix_type (str): Matrix type ('full' or 'jack')
            - title (str): Figure title (overriding the default)
        Returns:
            - fig: :mod:`matplotlib` figure object.
        
        Raises:
            Exception
        """
        if matrix_type=='full':
            if not hasattr(matrix_class,'reduced_matrix'):
                matrix_class.reduce_matrix('full')
            mat = matrix_class.reduced_matrix
            name = r'Reduced $C_{ab}$'
        elif matrix_type=='jack':
            if not hasattr(matrix_class,'reduced_jack_matrix'):
                matrix_class.reduce_matrix('jack')
            mat = matrix_class.reduced_jack_matrix
            name = r'Reduced $C^J_{ab}^J$'
        else:
            raise Exception("matrix_type must be 'full' or 'jack'")
        if len(title)>0:
            name = title
        return self._general_plotter(matrix_class.reduced_matrix,bin_min=bin_min,
                              bin_max=bin_max,title=name)
    def plot_precision(self,matrix_class,bin_min=None,bin_max=None):
        """ 
        Plot the precision matrix for specified bin ranges.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object
        
        Kwargs:
            - bin_min (int): Minimum bin number
            - bin_max (int): Maximum bin number
            - matrix_type (str): Matrix type ('full' or 'jack')
            - title (str): Figure title (overriding the default)
        Returns:
            - fig: :mod:`matplotlib` figure object.
        
        Raises:
            Exception
        """
        if matrix_type=='full':
            if not hasattr(matrix_class,'prec'):
                matrix_class.compute_precision('full')
            mat = matrix_class.prec
            name = r'$\Psi_{ab}$'
        elif matrix_type=='jack':
            if not hasattr(matrix_class,'prec_jack'):
                matrix_class.compute_precision('jack')
            mat = matrix_class.prec_jack
            name = r'$\Psi^J_{ab}^J$'
        else:
            raise Exception("matrix_type must be 'full' or 'jack'")
        if len(title)>0:
            name = title
        return self._general_plotter(matrix_class.prec,bin_min=bin_min,bin_max=bin_max,
                              title=name)
    def plot_reduced_precision(self,matrix_class,bin_min=None,bin_max=None):
        """ 
        Plot the precision matrix for specified bin ranges.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object
        
        Kwargs:
            - bin_min (int): Minimum bin number
            - bin_max (int): Maximum bin number
            - matrix_type (str): Matrix type ('full' or 'jack')
            - title (str): Figure title (overriding the default)
        Returns:
            - fig: :mod:`matplotlib` figure object.
        
        Raises:
            Exception
        """
        if matrix_type=='full':
            if not hasattr(matrix_class,'reduced_precision'):
                matrix_class.compute_precision('full')
            mat = matrix_class.reduced_precision
            name = r'Reduced $\Psi_{ab}$'
        elif matrix_type=='jack':
            if not hasattr(matrix_class,'reduced_precision_jack'):
                matrix_class.compute_precision('jack')
            mat = matrix_class.reduced_precision_jack
            name = r'Reduced $\Psi^J_{ab}^J$'
        else:
            raise Exception("matrix_type must be 'full' or 'jack'")
        if len(title)>0:
            name = title
        if not hasattr(matrix_class,'red_prec'):
            matrix_class.reduce_precision()
        return self._general_plotter(matrix_class.red_prec,bin_min,bin_max,
                              title=name)
    def plot_eigenvalues(self,matrix_class):
        """ Plot and compute the eigenvalues of the matrices. This computes the eigenvalues of both the full and jackknife matrices.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object 
            
        Returns:
            - figure: :mod:`matplotlib` figure
        """
        if not hasattr(matrix_class,'eigval'):
            matrix_class.compute_eigenvalues()
        figure=plt.figure()
        plt.plot(matrix_class.eigval,label='Full Matrix');
        plt.ylabel('Eigenvalue',fontsize=self.FS);plt.xlabel('Index',fontsize=self.FS)
        plt.plot(matrix_class.eigvalJ,label='Jackknife Matrix')
        plt.legend(fontsize=self.FS-4.);
        plt.yscale('log');
        return figure
    def plot_diagonal(self,matrix_class,fig=None,
                      use_jackknife=False,legend=False,name=""):
        """ Plot the diagonal elements of the full or jackknife matrices.
        
        Args:
            - matrix_class: :class:`CovarianceMatrix` object
        
        Kwargs:
            - fig: :mod:`matplotlib` figure object. The diagonal elements will be plotted onto this if provided, else a new figure will be created.
            - use_jackknife (bool): If True, plot the diagonal jackknife matrix elements. If False, plot the diagonal full matrix elements.
            - legend (bool): If True, include legend on plot.
            - name (str): String to include in legend (overwrites default)
            
        Returns:
            - fig: :mod:`matplotlib` figure object.
        """
        if fig==None:
            fig=plt.figure()
        ax=fig.gca()
        if use_jackknife:
            mat=matrix_class.c_jack_tot; ax.set_ylabel(r'$\mathrm{diag}(\hat{C}^J_{ab})$',fontsize=self.FS)
            if len(name)==0: 
                name='Jackknife Matrix'
        else:
            mat=matrix_class.c_tot; ax.set_ylabel(r'$\mathrm{diag}(\hat{C}_{ab})$',fontsize=self.FS)    
            if len(name)==0:
                name='Full Matrix'
        ax.plot(np.diag(mat),label=name); 
        ax.set_yscale('log'); ax.set_xlabel('Bin ID',fontsize=self.FS); ax.set_xlim([0,len(mat)])
        if legend:
            ax.legend(fontsize=self.FS-2)
        return fig
        
def KL_divergence(precision,covariance):
    """
    Utility function to return the negative log likelihood of the KL divergence between two matrices, one of which is inverted.
    The matrix known to higher precision should be inverted. 
    
    This has the following form: :math:`2\mathrm{KL}(\Psi,C) = \mathrm{Trace}(\Psi\times C) - \log\mathrm{det}\Psi - \log\mathrm{det}C - N_\mathrm{bins}`
    
    **NB**: The ordering of the precision and covariance matrices is non-trivial here.
    
    Args:
        - precision (np.ndarray): Precision matrix to be used
        - covariance (np.ndarray): Covariance matrix to be used
    
    Returns:
        - KL (float): Computed KL divergence
    
    """
    product = np.matmul(precision,covariance);
    N_bins=len(precision)
    logdetPrec = np.linalg.slogdet(precision)
    logdetCov = np.linalg.slogdet(covariance)
    if logdetPrec[0]!=1.:
        raise Exception('Undefined determinant for precision matrix')
    if logdetCov[0]!=1.:
        raise Exception('Undefined determinant for covariance matrix')
    KL = 0.5*(np.matrix.trace(product) - logdetPrec[1] - logdetCov[1] - N_bins)
    return KL
