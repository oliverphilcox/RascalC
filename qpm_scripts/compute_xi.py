## Convert coordinates of all QPM mocks
import subprocess
cwd = "/mnt/store1/oliverphilcox/DR12_QPM/"

for i in range(1,21):
    
    print("\nRUNNING FOR MOCK %d"%i)
    
    root= '/mnt/store1/oliverphilcox/DR12_QPM/'
    infile = root+'processed/qpm_galaxy_%s.xyzwj'%i

    subprocess.run(["python","/home/oliverphilcox/COMAJE/python/xi_estimator.py", "%s"%infile,"qpm_randoms_50x.xyzwj","qpm_randoms_10x.xyzwj","radial_binning_corr_fine.csv", "1.", "20", "20", "0", "xi_fine/mock_%s/"%i, "RRcounts_n90_m20.dat"],cwd=cwd)
    
    print("MOCK %d COMPLETE"%i)

print("All Mock Xi estimated")
