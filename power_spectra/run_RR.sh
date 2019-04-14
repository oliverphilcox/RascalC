#!/bin/bash

bash clean;
make;

for i in {16..50}
    do
        ./power -fname /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -fname2 /mnt/store1/oliverphilcox/PowerSpectra/qpm_randoms_50x_$i -out_string RR_$i
    done
