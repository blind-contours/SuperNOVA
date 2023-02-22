#!/bin/bash
# Job name:
#SBATCH --job-name=iter5_supernova
#
# Partition:
#SBATCH --partition=savio2
#

#SBATCH --qos=biostat_savio2_normal
#SBATCH --account=co_biostat

# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --exclusive
#
# Wall clock limit:
#SBATCH --time 72:00:00
#
## Command(s) to run (example):
module load r/4.0.3

### Run Simulation
R CMD BATCH --no-save ../03_run_simulation_5.R iter_5_supernova.Rout
