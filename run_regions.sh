#!/bin/bash

#PBS -l nodes=1:ppn=1

#PBS -l walltime=48:00:00

#PBS -A open

#PBS -m abe 	

#PBS -M tobias.gerken@psu.edu

# Get started

echo "Job started on ‘hostname‘ at ‘date‘"

# Load in matlab

module purge

module load matlab/R2019a

# Go to the correct place

cd $PBS_O_WORKDIR

# Run the job itself - a matlab script called runThis.m

matlab -nodisplay -nosplash -r "calculateBootstrap(1,1000)" > logRegion.matlabRun

# Finish up

echo "Job Ended at ‘date‘"