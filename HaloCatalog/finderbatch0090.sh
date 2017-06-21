#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalos
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH --mem=5G

#run the application:
srun -n 1 python ~/repositories/data-mining-galaxies/HaloCatalog/halo_finder0090.py
