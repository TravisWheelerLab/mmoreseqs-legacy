#!/usr/local/bin/bash
###########################################################################
#	NAME: 		mmseqs_plus.sh	
#	AUTHOR:		David Rich
#	DESC: 		Run main fb-pruner pipeline
###########################################################################

BENCH_DIR=$(pwd)/;
GENERAL_TMP=$BENCH_DIR/tmp/
FINAL_TMP=$GENERAL_TMP/$(uuidgen)/
TMP=$GENERAL_TMP/latest/

# check for valid command line arguments
if (( $# != 3 ))
then 
	echo "usage:  <target_file> <query_file> <dest_file>";
	exit;
fi

TARGET_FILE=$1;
QUERY_FILE=$2;
DEST_FILE=$3;

echo "target: $TARGET_FILE, query: $QUERY_FILE outfile: $DEST_FILE";

# build temporary if necessary
mkdir $GENERAL_TMP;
mkdir $FINAL_TMP;
mkdir $TMP;

# last step, move all temps