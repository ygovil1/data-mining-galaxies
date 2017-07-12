#!/bin/bash
#SBATCH -N 1
#SBATCH -J new_radiusfinderLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH --mem=4G

#run the application:
srun -n 1 python ./new_radiusfinder0070_500_1000.py
