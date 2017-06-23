#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG40110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 05:00:00
#SBATCH --mem=12G

#run the application:
srun -n 1 python ./halo_finder0110.py
