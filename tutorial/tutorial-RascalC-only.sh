#!/bin/bash
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128 # 128 hyperthreads = 64 physical cores
#SBATCH --job-name=RascalC-performance-test-tutorial

# load cosmodesi environment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

python -u tutorial-RascalC-only.py