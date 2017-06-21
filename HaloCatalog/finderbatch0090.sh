#!/bin/bash
#SBATCH -N 1
#SBATCH -J findhalos0090
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 05:50:00
#SBATCH --mem=16G

#run the application:
srun -n 1 python ~/repositories/data-mining-galaxies/HaloCatalog/halo_finder0090.py
