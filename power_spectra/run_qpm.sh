#!/bin/bash

bash clean;
make;

## Run DD counts
./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_galaxy_1.xyzwj -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_galaxy_1.xyzwj -out_string DD -output /mnt/store1/oliverphilcox/PowerQPM_New/R0_100/

## Run DR counts
./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_galaxy_1.xyzwj -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x.xyzwj -out_string DR -output /mnt/store1/oliverphilcox/PowerQPM_New/R0_100/

# 
# ## Run RR counts
# for i in {0..50}
#     do
#         ./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -out_string RR_$i -output /mnt/store1/oliverphilcox/PowerQPM50/
#     done
