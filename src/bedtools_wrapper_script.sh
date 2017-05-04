#!/bin/bash

## bedtools_wrapper_script.sh
## alex amlie-wolf 11-11-2015
## a wrapper script for bedtools to allow for calling bedtools from Python, while flexibly looking for
## a path to bedtools

BEDTOOLS_BIN_DIR=$1
shift
## the arguments are something like 'intersect -a ...'
## i.e. they do not include bedtools in the name
BEDTOOLS_ARGS=$*

## if we don't have bedtools, then try to find it
if [ "${BEDTOOLS_BIN_DIR}" == "None" ]; then
    ## first try just normal bedtools
    if ! type "bedtools" > /dev/null 2>&1; then
	## try loading bedtools or bedtools 2
	if type "module" > /dev/null 2>&1; then
	    module load bedtools > /dev/null 2>&1
	    ## if this still didn't work
	    if ! type "bedtools" > /dev/null 2>&1; then
		module load bedtools2 > /dev/null 2>&1
		if ! type "bedtools" > /dev/null 2>&1; then
		    echo "No bedtools found! Please insall bedtools and put it in your \$PATH, install it as a module, or provide a link to a bin/ directory for bedtools."
		    exit 1
		fi
	    fi
	fi
    fi
    ## if we reach this stage, then some bedtools should have worked    
    bedtools ${BEDTOOLS_ARGS}
else
    ${BEDTOOLS_BIN_DIR}/bedtools ${BEDTOOLS_ARGS}
fi
