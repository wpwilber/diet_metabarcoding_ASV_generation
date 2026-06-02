#!/bin/bash
#$ -M wwilber@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -q long
#$ -N merge_chim
#$ -cwd

source ~/.bash_profile
conda activate dada2-env

AMP=$1
Rscript scripts/merge_runs_remove_chimera.R "$AMP"

