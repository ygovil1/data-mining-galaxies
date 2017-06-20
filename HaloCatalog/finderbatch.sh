#!/bin/bash
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -J findhalos
#SBATCH --output=test.txt
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:01:00

#run the application:
srun -n 1  ~/repositories/data-mining-galaxies/HaloCatalog/prelim_halo_finder.py
