#!/bin/bash

#SBATCH --job-name="PBEC"
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH
srun ./PolaritonBEC
#srun ./venv/bin/python3 process.py
