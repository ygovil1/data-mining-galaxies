#!/bin/bash
#SBATCH -N 1
#SBATCH -J satellitefinderLG40090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 01:02:00
#SBATCH --mem=4G

#run the application:
srun -n 1 python ./satellite_finder.py
