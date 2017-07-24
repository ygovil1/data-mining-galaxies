#!/bin/bash
#SBATCH -N 1
#SBATCH -J satellitefinderLG38.0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 01:30:00
#SBATCH --mem=4G

#run the application:
srun -n 1 python ./satellitefinder.py
