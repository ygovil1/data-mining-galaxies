#!/bin/bash
#SBATCH -N 1
#SBATCH -J satellitefinderLG35.0070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00
#SBATCH --mem=4G

#run the application:
srun -n 1 python ./satellitefinder.py
