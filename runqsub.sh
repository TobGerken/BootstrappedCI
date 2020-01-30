#!/bin/bash

for CASES in 1 2 3; do
    for MODELS in 1 2; do
	for CATEGORIES in  1 2 3 4; do 
	       qsub -v CARG=$CASES,MARG=$MODELS,CATARG=$CATEGORIES runbatch.pbs 
	done
    done
done