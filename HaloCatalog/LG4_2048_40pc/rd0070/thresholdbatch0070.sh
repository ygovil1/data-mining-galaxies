#!/bin/bash
#SBATCH -N 1
#SBATCH -J thresh_finderLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 1:30:00
#SBATCH --mem=6G

#run the application:
srun -n 1 python ./threshold_finder.py
