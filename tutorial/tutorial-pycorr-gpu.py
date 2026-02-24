# Set the filenames
galaxies_filename = "mock_galaxy_DR12_CMASS_N_QPM_0001.rdzw"
randoms_filename = "mock_random_DR12_CMASS_N_50x1.rdzw"


# Read the original files: text with RA, DEC, Z (redshift) and weight columns.
import numpy as np
from astropy.table import Table
galaxies = Table(np.loadtxt(galaxies_filename, usecols = range(4)), names = ["RA", "DEC", "Z", "WEIGHT"]) # ignore the last column, not sure what it is
randoms = Table(np.loadtxt(randoms_filename, usecols = range(4)), names = ["RA", "DEC", "Z", "WEIGHT"])


# Compute the comoving distance within the fiducial (grid) cosmology. Here we use a utility function from the RascalC library to do this.
from RascalC.pre_process.convert_to_xyz import comoving_distance_Mpch
Omega_m = 0.29; Omega_k = 0; w_DE = -1 # density parameters of matter and curvature, and the equation-of-state parameter of dark energy
galaxies["comov_dist"] = comoving_distance_Mpch(galaxies["Z"], Omega_m, Omega_k, w_DE)
randoms["comov_dist"] = comoving_distance_Mpch(randoms["Z"], Omega_m, Omega_k, w_DE)


# Let us define a utility function for position formatting that will be useful on several more occasions.
def get_rdd_positions(catalog: Table) -> tuple[np.typing.NDArray[np.float64]]: # utility function to format positions from a catalog
    return (catalog["RA"], catalog["DEC"], catalog["comov_dist"])


# Assign jackknife regions to both galaxies and randoms.
from RascalC.pre_process.create_jackknives_pycorr import get_subsampler_xirunpc
n_jack = 60 # number of regions
subsampler = get_subsampler_xirunpc(get_rdd_positions(galaxies), n_jack, position_type = "rdd") # "rdd" means RA, DEC in degrees and then distance (corresponding to pycorr)
galaxies["JACK"] = subsampler.label(get_rdd_positions(galaxies), position_type = "rdd")
randoms["JACK"] = subsampler.label(get_rdd_positions(randoms), position_type = "rdd")


# Select a smaller subset of randoms to make pair counting and `RascalC` importance sampling more feasible
x_randoms = 10 # how many times the number of galaxies should the number of randoms be; the total number of randoms is ≈50x the number of galaxies
np.random.seed(42) # for reproducibility
randoms_subset = randoms[np.random.choice(len(randoms), x_randoms * len(galaxies), replace = False, p = randoms["WEIGHT"] / np.sum(randoms["WEIGHT"]))]


# ## Pair counts and correlation functions with [`pycorr`](https://github.com/cosmodesi/pycorr)

# First, we choose whether to use full randoms or a smaller subset for pair counting.
# The latter is faster but a bit less precise.
randoms_for_counts = randoms
# randoms_for_counts = randoms_subset


# We continue by splitting the randoms into parts of roughly the same size as data.
# This gives high precision at fixed computing cost [(Keihänen et al 2019)](https://arxiv.org/abs/1905.01133).
n_splits = int(np.rint(len(randoms_for_counts) / len(galaxies))) # the number of parts to split the randoms to
print(f"Splitting randoms into {n_splits} parts")

# split randoms into the desired number of parts randomly
random_indices = np.arange(len(randoms_for_counts))
np.random.seed(42) # for reproducibility
np.random.shuffle(random_indices) # random shuffle in place
random_parts = [randoms_for_counts[random_indices[i_random::n_splits]] for i_random in range(n_splits)] # quick way to produce parts of almost the same size

# normalize the weights in each part — fluctuations in their sums may be a bit of a problem
for i_random in range(n_splits): random_parts[i_random]["WEIGHT"] /= np.sum(random_parts[i_random]["WEIGHT"])


# Import libraries and select settings
from pycorr import TwoPointCorrelationFunction, setup_logging
from tqdm import trange # nice progress bar
import os


n_threads = 4 # 4 GPUs
gpu = True
split_above = 20 # Mpc/h. Below this, will use concatenated randoms. Above, will use split.
s_max = 200 # maximal separation in Mpc/h
n_mu = 200 # number of angular (µ) bins
counts_filename = f"allcounts_mock_galaxy_DR12_CMASS_N_QPM_0001_lin_njack{n_jack}_nran{n_splits}_split{split_above}.npy" # filename to save counts


mu_edges = np.linspace(-1, 1, n_mu + 1) # make uniform µ bins between -1 and 1, or twice less bins between 0 and 1 after wrapping (will be done within RascalC wrapper)
s_edges_all = (np.arange(split_above + 1), np.arange(split_above, s_max + 1)) # 1 Mpc/h wide separation bins from 0 to s_max Mpc/h, separated to concatenated/split random regions. Can be rebinned to any bin width that divides split_above and s_max

def run_pair_counts(): # the code to run pycorr, needs to be invoked in a separate process from what will run RascalC then!
    setup_logging()
    results = []
    # compute
    for i_split_randoms, s_edges in enumerate(s_edges_all):
        result = 0
        D1D2 = None # to compute the data-data counts on the first go but not recompute then
        for i_random in trange(n_splits if i_split_randoms else 1, desc="Computing counts with random part"):
            these_randoms = random_parts[i_random] if i_split_randoms else randoms_for_counts
            tmp = TwoPointCorrelationFunction(mode = 'smu', edges = (s_edges, mu_edges),
                                            data_positions1 = get_rdd_positions(galaxies), data_weights1 = galaxies["WEIGHT"], data_samples1 = galaxies["JACK"],
                                            randoms_positions1 = get_rdd_positions(these_randoms), randoms_weights1 = these_randoms["WEIGHT"], randoms_samples1 = these_randoms["JACK"],
                                            position_type = "rdd", engine = "corrfunc", D1D2 = D1D2, gpu = gpu, nthreads = n_threads)
            # "rdd" means RA, DEC in degrees and then distance
            D1D2 = tmp.D1D2 # once computed, becomes not None and will not be recomputed
            result += tmp
        results.append(result)
    corr = results[0].concatenate_x(*results) # join the unsplit and split parts
    corr.D1D2.attrs['nsplits'] = n_splits

    corr.save(counts_filename)


# run pair counts
run_pair_counts()
