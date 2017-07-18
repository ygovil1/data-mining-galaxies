#!/bin/bash
#SBATCH -N 1
#SBATCH -J ratiofinderLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./ratio_finder0070.py
