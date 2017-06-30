#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG76_0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 05:00:00
#SBATCH --mem=10G

#run the application:
srun -n 1 python ./halo_finder0110.py
