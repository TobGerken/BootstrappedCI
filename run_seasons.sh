#!/bin/bash

#PBS -l nodes=1:ppn=1

#PBS -l walltime=5:00

#PBS -A open

# Get started

echo "Job started on ‘hostname‘ at ‘date‘"

# Load in matlab

module purge

module load matlab/R2016a

# Go to the correct place

cd $PBS_O_WORKDIR

# Run the job itself - a matlab script called runThis.m

matlab-bin -nodisplay -nosplash < calculateBootstrap({'Seasons'},1000) > logSeasons.matlabRun

# Finish up

echo "Job Ended at ‘date‘"