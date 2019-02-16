## Convert coordinates of all QPM mocks
import subprocess
cwd = "/home/oliverphilcox/COMAJE/"

for i in range(1,21):
    
    print("\nRUNNING FOR MOCK %d"%i)
    
    root= '/mnt/store1/oliverphilcox/DR12_QPM/'
    infile = root+'mock_galaxy_DR12_CMASS_N_QPM_%s.rdzw'%(str(i).zfill(4))
    outfile = root+'processed/qpm_galaxy_%s.xyz'%i
    
    subprocess.run(["python","python/convert_to_xyz.py", "%s"%infile,"%s"%outfile,"0.29", "0.", "-1."],cwd=cwd)
    
    outfile2 = root+'processed/qpm_galaxy_%s.xyzwj'%i

    print("Creating Jackknives for mock %d"%i)
    subprocess.run(["python","python/create_jackknives.py","%s"%outfile,"%s"%outfile2,"8"],cwd=cwd)
    
    print("MOCK %d COMPLETE"%i)
