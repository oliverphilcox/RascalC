## Convert coordinates of all QPM mocks
import subprocess
cwd = "/mnt/store1/oliverphilcox/DR12_QPM/"

for i in range(1,21):
    
    print("\nRUNNING FOR MOCK %d"%i)
    
    root= '/mnt/store1/oliverphilcox/DR12_QPM/'
    infile = root+'processed/qpm_galaxy_%s.xyzwj'%i

    subprocess.run(["python","/home/oliverphilcox/COMAJE/python/xi_estimator_jack.py", "%s"%infile,"qpm_randoms_50x.xyzwj","qpm_randoms_10x.xyzwj","radial_binning_cov.csv", "1.", "10", "20", "0", "xi_jack/mock_%s/"%i, "jackknife_pair_counts_n35_m10_j169_11.dat"],cwd=cwd)
    
    print("MOCK %d COMPLETE"%i)

print("All Mock Xi-Jack estimated")
