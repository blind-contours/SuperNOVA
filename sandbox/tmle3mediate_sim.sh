#!/bin/bash
# Job name:
#SBATCH --job-name=tmle3mediate-simulation
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run (example):
module load r/3.6.3
module load r-packages

## Set SBATCH_ACCOUNT and SBATCH_QOS env variables prior to running

### for foreach+doSNOW ###
R CMD BATCH --no-save 03_run_simulation.R
