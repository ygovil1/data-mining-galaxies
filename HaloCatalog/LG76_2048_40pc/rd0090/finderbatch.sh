#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG76_0090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 10:00:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./halo_finder.py
