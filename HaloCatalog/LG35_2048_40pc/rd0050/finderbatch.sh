#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalosLG35.0050
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 05:30:00
#SBATCH --mem=10G

#run the applications
srun -n 1 python ./halo_finder.py
srun -n 1 python ./list_maker.py

