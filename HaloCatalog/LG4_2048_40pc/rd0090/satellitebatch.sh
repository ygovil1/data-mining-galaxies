#!/bin/bash
#SBATCH -N 1
#SBATCH -J satellitefinderv2LG40090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:59:00
#SBATCH --mem=4G

#run the application:
srun -n 1 python ./satellite_finder.py