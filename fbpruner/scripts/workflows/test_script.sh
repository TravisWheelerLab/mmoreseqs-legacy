#!/usr/bin/bash
###########################################################################
#	NAME: 		test_script.sh	
#	AUTHOR:		David Rich
#	DESC: 		runs mmseqs, then pipes output to fb-pruner search
###########################################################################

# SET VARIABLES

# commandline variables (verify proper number of variables)
NUM_ARGS=$#
if (( $NUM_ARGS < 2 )); then 
	echo "ERROR: inncorrect number of parameters: $NUM_ARGS of 2"
	echo "Usage: <TARGET> <QUERY>"
	exit
fi

# -- Main Args -- #
TARGET=$1
QUERY=$2

echo "This is a test_script.sh"
echo "=== ARGS ==="
echo "TARGET: $TARGET"
echo "QUERY: $QUERY"
