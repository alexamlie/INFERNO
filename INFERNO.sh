#!/bin/bash

## INFERNO.sh
## Alex Amlie-Wolf
## wrapper script for the main INFERNO python function so that you can run it in bsub

## for now, hard-code and load the appropriate modules
## TODO: let users define these to make sure everything is loaded
module load python/2.7.9
module load bedtools2
module load R/3.2.3
module load plink/1.90Beta

## just feed all the arguments to the python script
python ./INFERNO.py $*
