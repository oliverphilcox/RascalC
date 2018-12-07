def read_RR_all(n,m,file_root,a,string='full'):
    """
    Read in the estimated and true RR pair counts.
    n = no. radial bins
    m = no. angular bins
    file_root = directory housing files
    a = shot-noise rescaling parameter
    string = defines which estimate to use. 'full' gives the total estimate
    """
    rrfile=file_root+'RR_n%d_m%d_%s.txt' %(n,m,string)
    rr_est = np.loadtxt(rrfile)*corr_factor**2.
    rr_true = np.loadtxt(rr_true_file)
    return rr_est, rr_true

def full_covariance(n,m,file_root,a,string='full'):
    """
    Read in the full covariance matrix estimates.
    n = no. radial bins
    m = no. angular bins
    file_root = directory housing files
    a = shot-noise rescaling parameter
    string = defines which estimate to use. 'full' gives the total estimate
    """
    c2file=file_root+'c2_n%d_m%d_%s.txt' %(n,m,string)
    c3file=file_root+'c3_n%d_m%d_%s.txt' %(n,m,string)
    c4file=file_root+'c4_n%d_m%d_%s.txt' %(n,m,string)
    errc4file=file_root+'c4err_n%d_m%d_%s.txt' %(n,m,string)
    c2=np.diag(np.loadtxt(c2file))*corr_factor**2.
    c3=np.loadtxt(c3file)*corr_factor**3.
    c4=np.loadtxt(c4file)*corr_factor**4.
    errc4=np.loadtxt(errc4file)*corr_factor**8.
    c3=(c3+c3.T)/2.
    c4=(c4+c4.T)/2.
    errc4=(errc4+errc4.T)/2.
    c_tot=c4+a**2.*c2+a*c3
    return c2,c3,c4,errc4,c_tot
    
def jackknife_covariance(n,m,file_root,a,string='full'):
    """
    Read in the jackknife covariance matrix estimates.
    n = no. radial bins
    m = no. angular bins
    file_root = directory housing files
    a = shot-noise rescaling parameter
    string = defines which estimate to use. 'full' gives the total estimate
    """
    c2file=file_root+'c2j_n%d_m%d_%s.txt' %(n,m,string)
    c3file=file_root+'c3j_n%d_m%d_%s.txt' %(n,m,string)
    c4file=file_root+'c4j_n%d_m%d_%s.txt' %(n,m,string)
    cxfile=file_root+'cxj_n%d_m%d_%s.txt' %(n,m,string)
    errc4file=file_root+'c4errj_n%d_m%d_%s.txt' %(n,m,string)
    errcxfile=file_root+'cxerrj_n%d_m%d_%s.txt' %(n,m,string)
    c2=np.diag(np.loadtxt(c2file))*corr_factor**2.
    c3=np.loadtxt(c3file)*corr_factor**3.
    c4=np.loadtxt(c4file)*corr_factor**4.
    cx=np.loadtxt(cxfile)*corr_factor**4.
    errc4=np.loadtxt(errc4file)*corr_factor**8.
    c3=(c3+c3.T)/2.
    c4=(c4+c4.T)/2.
    cx=(cx+cx.T)/2.
    errc4=(errc4+errc4.T)/2.
    c_tot=c4+a**2.*c2+a*c3+cx
    return c2,c3,c4,cx,errc4,c_tot
