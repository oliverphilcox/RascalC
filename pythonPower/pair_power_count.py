import numpy as np
from scipy.special import spherical_jn

### PARAMETERS

gal_name = '/mnt/store1/oliverphilcox/Legendre2PCF/qpm_galaxy_1.xyzwj'
R0=200 # max radius in Mpc/h
Rcut=100 # transition radius for window function
k_cen = np.arange(0.1,1.1,0.1)

# Load particles
print("Loading particles")
gal_file=np.loadtxt(gal_name)

# Initialize things
gal_x=gal_file[:,0]
gal_y=gal_file[:,1]
gal_z=gal_file[:,2]
gal_w=gal_file[:,3]
N_gal=len(gal_x)
R02=np.power(R0,2.)
k_vec = np.reshape(k_cen,(1,-1))
b=Rcut/R0
a=(R0-Rcut)/(2.*R0)
pair_counts = np.zeros(len(k_cen))
print("Loaded %d particles"%N_gal)

## Define window function
def windower(r):
    x = r/R0
    output = np.zeros_like(r)
    filt1 = np.where(x<b)
    filt2 = np.where((x>b)&(x<b+a))
    filt3 = np.where((x>b+a)&(x<1))
    output[filt1] = 1.
    output[filt2] = 1.-((x[filt2]-b)/a)**3.+0.5*((x[filt2]-b)/a)**4.
    output[filt3] = -((x[filt3]-b-2*a)/a)**3.-0.5*((x[filt3]-b-2*a)/a)**4.
    return output

# Create runner
def run_1000(args):
    start,stop=args
    pair_counts=np.zeros(len(k_cen))
    for i in range(start,stop):
        this_x=gal_x[i]
        this_y=gal_y[i]
        this_z=gal_z[i]
        this_w=gal_w[i]
        filt1=np.where((np.abs(gal_x-this_x)<R0)&(np.abs(gal_y-this_y)<R0)&(np.abs(gal_z-this_z)<R0))
        rel_dis2 = np.power(gal_x[filt1]-this_x,2.)+np.power(gal_y[filt1]-this_y,2.)+np.power(gal_z[filt1]-this_z,2.)
        filt2=np.where(rel_dis2<R02)
        # Select good galaxies
        good_x=gal_x[filt1][filt2]
        good_y=gal_y[filt1][filt2]
        good_z=gal_z[filt1][filt2]
        rel_dis = np.sqrt(rel_dis2[filt2])
        k_dis = np.matmul(np.reshape(rel_dis,(-1,1)),k_vec)
        weight_prod = np.reshape(gal_w[filt1][filt2]*windower(rel_dis),(-1,1))
        conts = np.multiply(spherical_jn(0,k_dis),weight_prod)
        pair_counts+=np.sum(conts,axis=0)
    return pair_counts

if __name__=='__main__':
    thousand_list=[]
    for i in range(N_gal):
        if i%1000==0:
            thousand_list.append([i,min(i+1000,N_gal)])
    import multiprocessing as mp
    p=mp.Pool(20)
    import tqdm
    
    print("Starting multiprocessing")
    output=list(tqdm.tqdm(p.imap_unordered(run_1000,thousand_list),total=len(thousand_list)))
    print("Multiprocessing complete")
    
    full_pair_counts = np.zeros(len(k_cen))
    for o in output:
        full_pair_counts+=o
    
    np.savez('/mnt/store1/oliverphilcox/PowerSpectra/binned_k_power.npz',k=k_cen,counts=full_pair_counts)
