#!/bin/bash
#SBATCH --account=desi_g
#SBATCH --constraint=gpu
#SBATCH --qos=regular
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --job-name=RascalC-tutorial-pycorr

# load cosmodesi environment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

python -u tutorial-pycorr-gpu.py