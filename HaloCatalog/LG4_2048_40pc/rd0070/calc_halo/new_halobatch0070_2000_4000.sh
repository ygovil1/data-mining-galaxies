#!/bin/bash
#SBATCH -N 1
#SBATCH -J new_radiusfinderLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 05:59:00
#SBATCH --mem=8G

#run the application:
srun -n 1 python ./new_radiusfinder0070_2000_4000.py
