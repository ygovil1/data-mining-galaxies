#!/bin/bash
#SBATCH -N 1
#SBATCH -J agefinderLG38.0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 01:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./agefinder.py
