#!/bin/bash
#SBATCH -N 1
#SBATCH -J new_radiusfinderLG35.0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 10:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./new_radiusfinder.py
