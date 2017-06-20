#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalos
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 03:00:00

#run the application:
srun -n 1 python ~/repositories/data-mining-galaxies/HaloCatalog/prelim_halo_finder.py
