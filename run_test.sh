#!/bin/bash

#PBS -l nodes=1:ppn=1

#PBS -l walltime=20:00

#PBS -A open

# Get started

echo "Job started on ‘hostname‘ at ‘date‘"

# Load in matlab

module purge

module load matlab/R2019a

# Go to the correct place

cd $PBS_O_WORKDIR

# Run the job itself - a matlab script called runThis.m

matlab -nodisplay -nosplash -r "CalculateBootstrap $CARG 1 1 1 1" > logTest.matlabRun

# Finish up

echo "Job Ended at ‘date‘"
