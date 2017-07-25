#!/bin/bash
#SBATCH -N 1
#SBATCH -J stelradfinderLG4.0070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 12:59:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./stel_radfinder0070.py
