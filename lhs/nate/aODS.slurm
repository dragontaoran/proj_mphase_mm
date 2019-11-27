#!/bin/bash

# Beginning of SBATCH batch script.
#SBATCH --mail-user=jonathan.schildcrout@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output=myjob.output
#SBATCH --array=1-10
##SBATCH --job-name=lhs_example

HOMEDIR=~/rsch/OutDepSamp/Projects/AODSBinary/LHSAnalysis/ACCRE

mkdir -p $HOMEDIR/results/seed_$SLURM_ARRAY_TASK_ID
mkdir -p $HOMEDIR/tmp/seed_$SLURM_ARRAY_TASK_ID
mkdir -p $HOMEDIR/lhs/seed_$SLURM_ARRAY_TASK_ID
mkdir -p $HOMEDIR/adaptive_ss/seed_$SLURM_ARRAY_TASK_ID
mkdir -p $HOMEDIR/adaptive_des/seed_$SLURM_ARRAY_TASK_ID
mkdir -p $HOMEDIR/fixed2stage/seed_$SLURM_ARRAY_TASK_ID


# substitute parameters in R file and save it
sed s/bbb/$SLURM_ARRAY_TASK_ID/g $HOMEDIR/aODS_April2018.R >$HOMEDIR/tmp/seed_$SLURM_ARRAY_TASK_ID/lhs_$SLURM_ARRAY_TASK_ID.R

# Now submit the slurm script
cd $HOMEDIR/tmp/seed_$SLURM_ARRAY_TASK_ID/

module purge
module load GCC/5.4.0-2.26
module load OpenMPI/1.10.3
module load R/3.3.3-X11-20160819

R CMD BATCH $HOMEDIR/tmp/seed_$SLURM_ARRAY_TASK_ID/lhs_$SLURM_ARRAY_TASK_ID.R  $HOMEDIR/tmp/seed_$SLURM_ARRAY_TASK_ID/lhs_$SLURM_ARRAY_TASK_ID.Rout
