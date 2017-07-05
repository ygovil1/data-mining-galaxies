#!/bin/bash
#SBATCH -N 1
#SBATCH -J rad_ratiofinderLG40070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 4:30:00
#SBATCH --mem=20G

#run the application:
srun -n 1 python ./rad_ratiofinder0070.py
