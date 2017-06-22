#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalos0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 30:00:00
#SBATCH --mem=60G

#run the application:
srun -n 1 python ./RD0110/halo_finder0110.py
