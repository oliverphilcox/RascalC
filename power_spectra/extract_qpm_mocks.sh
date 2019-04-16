#!/bin/bash

for i in {0001..0100}
    do
        echo "Converting file $i"
        python ~/COMAJE/python/convert_to_xyz.py /mnt/store1/oliverphilcox/DR12_QPM/unprocessed/mock_galaxy_DR12_CMASS_N_QPM_$i.rdzw /mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_$i.xyzw 0.29 0. -1. 
    done
