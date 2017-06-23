#!/bin/bash
#SBATCH -N 1
#SBATCH -J makeplots0070
#SBATCH --mail-user=ygovil@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:30:00
#SBATCH --mem=8G

#run the application:
srun -n 1 python ./make_plot0070.py
