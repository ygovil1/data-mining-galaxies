#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 5:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./halo_finder0070.py
