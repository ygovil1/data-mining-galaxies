#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 16:00:00
#SBATCH --mem=60G

#run the application:
srun -n 1 python ./ratio_finder0070.py
