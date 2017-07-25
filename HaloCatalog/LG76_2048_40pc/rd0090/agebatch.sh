#!/bin/bash
#SBATCH -N 1
#SBATCH -J agefinderLG76.0090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 02:29:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./agefinder.py
