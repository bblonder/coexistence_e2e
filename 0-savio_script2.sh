#!/bin/bash
# Email me:
#SBATCH --mail-user=benjamin.blonder@berkeley.edu
# Job name:
#SBATCH --job-name=ben_love
#
# Partition:
#SBATCH --partition=savio4_htc
#
# QoS:
#SBATCH --qos=savio_normal
#
# Account:
#SBATCH --account=fc_mel
#
# Request one node:
#SBATCH --nodes=1
#
# Number of processors for threading:
#SBATCH --cpus-per-task=20
#
# Wall clock limit:
#SBATCH --time=36:00:00
#
## Command(s) to run:
module load r
module load r-packages
module load r-spatial
cd /global/scratch/users/benjaminblonder/love_august
R CMD BATCH --no-save 1-predict-test.R job.Rout