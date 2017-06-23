#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG40050
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 04:30:00
#SBATCH --mem=10G

#run the application:
srun -n 1 python ./halo_finder0050.py
