### Function to fit a model to the survey correction function, defined as the ratio between model and true RR pair counts for a single survey. This fits a piecewise polynomial model to the data.

## NB: Input RR counts should NOT be normalized by summed galaxy weights here.
## NB: Assume mu is in [0,1] limit here

import os
import numpy as np
import scipy.spatial as ss
from scipy.optimize import curve_fit


def compute_phi_periodic(w_bar1: float, w_bar2: float, n_bar1: float, n_bar2: float, V: float, n: int) -> np.ndarray[float]:
    ## Define normalization constant
    this_norm = n_bar1 * n_bar2 * w_bar1 * w_bar2 # keep it simple, only need to be self-consistent - use the same norm everywhere

    norm = V * this_norm

    phi = np.zeros((n, 7))
    phi[:, 0] = 1.
    phi[:, 3] = 1.

    return phi / norm


def compute_phi_aperiodic(w_bar1: float, w_bar2: float, n_bar1: float, n_bar2: float, V: float, r_bins: np.ndarray[float], this_RR: np.ndarray[float], index: str, print_function = print) -> np.ndarray[float]:
    print_function("\nComputing survey correction factor %s"%index)
    ## Define normalization constant
    this_norm = n_bar1 * n_bar2 * w_bar1 * w_bar2 # keep it simple, only need to be self-consistent - use the same norm everywhere

    norm = V * this_norm

    # Find radial bin volumes
    vol_r = 4.*np.pi/3.*(r_bins[:,1]**3 - r_bins[:,0]**3)

    n, m = this_RR.shape

    # Generate mu bins
    mu_cen = np.arange(1./(2.*m),1.+1./(2.*m),1./m)
    delta_mu = mu_cen[-1]-mu_cen[-2]
    assert(m==len(mu_cen))

    ## Define RR model
    def RR_model(r_bin, mu):
        return norm * vol_r[r_bin] * delta_mu

    # Compute correction functions
    Phi_values = []
    for r_bin in range(n):
        Phi_values.append(RR_model(r_bin, mu_cen)/this_RR[r_bin,:])

    ## check for order of magnitude consistency
    if np.mean(Phi_values) > 1e2:
        raise Exception("RR_true seems to be much smaller than RR_model. Is the input RR consistent with weights (not normalized)?")
    if np.mean(Phi_values) < 1e-2:
        raise Exception("RR_true seems to be much larger than RR_model. Is the input RR consistent with weights (not normalized)?")

    ## Define Phi model (piecewise continuous polynomial)
    mu_crit=0.75
    def Phi_model(mu,a0,a1,a2,b2,b3,mu_crit=mu_crit):
        b1=a1+2*mu_crit*(a2-b2)-3*b3*mu_crit**2.
        b0=a0+(a1-b1)*mu_crit+(a2-b2)*mu_crit**2.-b3*mu_crit**3.
        filt1=np.where(mu<mu_crit)
        filt2=np.where(mu>=mu_crit)
        output=np.zeros_like(mu)
        output[filt1]=a0+a1*mu[filt1]+a2*mu[filt1]**2.
        output[filt2]=b0+b1*mu[filt2]+b2*mu[filt2]**2.+b3*mu[filt2]**3.
        return output

    ## Find optimal parameters
    def fit_model(mu,good_param):
        return Phi_model(mu,good_param[0],good_param[1],good_param[2],good_param[3],good_param[4])

    fit_params=[]
    errors=[]
    for i in range(n):
        good_param,_=curve_fit(Phi_model,mu_cen,Phi_values[i],p0=[0,0,0,0,0])
        a0,a1,a2,b2,b3=good_param
        b1=a1+2*mu_crit*(a2-b2)-3*b3*mu_crit**2.
        b0=a0+(a1-b1)*mu_crit+(a2-b2)*mu_crit**2.-b3*mu_crit**3.
        out_params=[a0,a1,a2,b0,b1,b2,b3]
        fit_params.append(out_params)
        errors.append(np.abs(Phi_values[i]-fit_model(mu_cen,good_param))/fit_model(mu_cen,good_param))

    print_function("Fitted dataset %s with mean fractional error %.1e"%(index,np.mean(errors)))

    return np.asarray(fit_params) / norm

def compute_V_n_w_bar(randoms_pos: np.ndarray[float], gal_w: np.ndarray[float]) -> tuple[float, float, float]:
    w_bar = np.mean(gal_w)
    N_gal = len(gal_w)

    if randoms_pos.shape[1] != 3: randoms_pos = randoms_pos.T
    if randoms_pos.shape[1] != 3: raise ValueError("Galaxy positions do not appear 3D")

    ## Find survey volume via ConvexHull in Scipy
    V = ss.ConvexHull(randoms_pos).volume # presumably in (Mpc/h)^3

    n_bar = N_gal / V
    return V, n_bar, w_bar

def compute_V_n_w_bar_from_file(random_file: str, index = 1, print_function = print) -> tuple[float, float, float]:
    print_function(f"Loading galaxy set {index}")
    all_gal = np.loadtxt(random_file)

    V, n_bar, w_bar = compute_V_n_w_bar(all_gal[:, :3], all_gal[:, 3])
    print_function(f'Survey volume {index} is approximately: {V/1e9:.2f} (Gpc/h)^3')
    return V, n_bar, w_bar


def load_RR(RR_file: str, n: int) -> np.ndarray[float]:
    RR_flat = np.loadtxt(RR_file) # not change normalization here
    return RR_flat.reshape((n, -1))


def compute_correction_function(random_file: str, binfile: str, outdir: str, periodic: bool, RR_file: str | None = None, print_function = print) -> None:
    if periodic:
        print("Assuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere")
    elif RR_file is None: raise TypeError("The RR file must be specified if aperiodic")
    
    V, n_bar, w_bar = compute_V_n_w_bar(random_file)

    # Load in binning files 
    r_bins = np.loadtxt(binfile)
    n = len(r_bins)

    if periodic:
        ## Output periodic survey correction function
        phi = compute_phi_periodic(w_bar, w_bar, n_bar, n_bar, V, n)

    else:
        ## load in RR counts
        RR = load_RR(RR_file, n)
        m = RR.shape[1]

        phi = compute_phi_aperiodic(w_bar, w_bar, n_bar, n_bar, V, r_bins, RR, print_function = print_function)

    if periodic:
        outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_periodic_11.txt'%(n))
    else:
        outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_m%d_11.txt'%(n, m))
        
    np.savetxt(outfile, phi)
    print_function("Saved (normalized) output to %s\n" % outfile)


def compute_correction_function_multi(random_file: str, random_file2: str, binfile: str, outdir: str, periodic: bool, RR_file: str | None = None, RR_file12: str | None = None, RR_file2: str | None = None, print_function = print) -> None:
    if periodic:
        print_function("Assuming periodic boundary conditions - so Phi(r,mu) = 1 everywhere")
    elif any(file is None for file in (RR_file, RR_file12, RR_file2)):
        raise TypeError("All the RR files must be specified if aperiodic")

    V1, n_bar1, w_bar1 = compute_V_n_w_bar_from_file(random_file, print_function = print_function)
    V2, n_bar2, w_bar2 = compute_V_n_w_bar_from_file(random_file2, index = 2, print_function = print_function)

    # Load in binning files
    r_bins = np.loadtxt(binfile)
    n = len(r_bins)

    if periodic:
        ## Periodic correction function is simple

        phi_11 = compute_phi_periodic(w_bar1, w_bar1, n_bar1, n_bar1, V1, n)
        phi_12 = compute_phi_periodic(w_bar1, w_bar2, n_bar1, n_bar2, np.sqrt(V1*V2), n)
        phi_22 = compute_phi_periodic(w_bar2, w_bar2, n_bar2, n_bar2, V2, n)
    else:
        ## Continue for aperiodic case

        # Load RR counts
        RR = load_RR(RR_file, n)
        m = RR.shape[1]

        RR2 = load_RR(RR_file2, n)
        if RR2.shape != RR.shape: raise ValueError("Need the same binning for all RR files")

        RR12 = load_RR(RR_file12, n)
        if RR12.shape != RR.shape: raise ValueError("Need the same binning for all RR files")

        # Now run this
        phi_11 = compute_phi_aperiodic(w_bar1, w_bar1, n_bar1, n_bar1, V1, r_bins, RR, "11", print_function)
        phi_12 = compute_phi_aperiodic(w_bar1, w_bar2, n_bar1, n_bar2, np.sqrt(V1*V2), r_bins, RR, "12", print_function) # use geometrically average volume, should be not critical
        phi_12 = compute_phi_aperiodic(w_bar2, w_bar2, n_bar2, n_bar2, V2, RR, r_bins, "22", print_function)

    roots = ['11', '12', '22']
    phis = [phi_11, phi_12, phi_22]

    for phi, index in zip(phis, roots):
        if periodic:
            outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_periodic_%s.txt'%(n, index))
        else:
            outfile = os.path.join(outdir, 'BinCorrectionFactor_n%d_m%d_11.%s'%(n, m, index))
        np.savetxt(outfile, phi)
        print_function("Saved (normalized) output for field %s to %s"%(index, outfile))