#!/bin/bash

bash clean;
make;

#     
# ## Run RR counts
# for i in {0..50}
#     do
#         ./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -out_string RR_$i -output /mnt/store1/oliverphilcox/PowerQPM/
#     done

for i in {0001..0030}
    do 
        ## Run DD counts
        ./power -fname /mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_$i.xyzw -fname2 /mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_$i.xyzw -out_string DD_mock_$i -output /mnt/store1/oliverphilcox/PowerQPM_New/R0_100_noPhi/
        
        ## Run DR counts
        ./power -fname /mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_$i.xyzw -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x.xyzwj -out_string DR_mock_$i -output /mnt/store1/oliverphilcox/PowerQPM_New/R0_100_noPhi/

    done

