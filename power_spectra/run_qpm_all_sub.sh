#!/bin/bash

bash clean;
make;

#     
# ## Run RR counts
# for i in {0..50}
#     do
#         ./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -out_string RR_$i -output /mnt/store1/oliverphilcox/PowerQPM/
#     done

for i in {0050..0100}
    do 
        ## Run DD counts
        ./power -fname /mnt/store1/oliverphilcox/SubsampledPower/qpm_galaxy_${i}_sub.xyzw -fname2 /mnt/store1/oliverphilcox/SubsampledPower/qpm_galaxy_${i}_sub.xyzw -out_string DD_mock_$i -R0 200 -nside 51 -binfile /mnt/store1/oliverphilcox/SubsampledPower/k_binning_sub.csv -output /mnt/store1/oliverphilcox/SubsampledPower/R0_200/
    
         ## Run DR counts
        ./power -fname /mnt/store1/oliverphilcox/SubsampledPower/qpm_galaxy_${i}_sub.xyzw -fname2 /mnt/store1/oliverphilcox/SubsampledPower/qpm_randoms_50x_sub.xyzwj -out_string DR_mock_$i -R0 200 -nside 101 -binfile /mnt/store1/oliverphilcox/SubsampledPower/k_binning_sub.csv -output /mnt/store1/oliverphilcox/SubsampledPower/R0_200/

    done

