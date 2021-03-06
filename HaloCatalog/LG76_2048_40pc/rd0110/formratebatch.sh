#!/bin/bash
#SBATCH -N 1
#SBATCH -J formrateLG76.0110
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 2:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./formratefinder.py
