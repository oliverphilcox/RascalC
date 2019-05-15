#!/bin/bash

bash clean;
make;

for i in {0..50}
    do
        ./power -fname /mnt/store1/oliverphilcox/SubsampledPower/qpm_randoms_50x_${i}_sub -fname2 /mnt/store1/oliverphilcox/SubsampledPower/qpm_randoms_50x_${i}_sub -out_string RR_$i -R0 200 -nside 51 -binfile /mnt/store1/oliverphilcox/SubsampledPower/k_binning_sub.csv -output /mnt/store1/oliverphilcox/SubsampledPower/R0_200/
    done
