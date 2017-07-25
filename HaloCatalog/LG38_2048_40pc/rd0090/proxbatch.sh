#!/bin/bash
#SBATCH -N 1
#SBATCH -J proxfinderLG38.0090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 03:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./proxfinder.py
