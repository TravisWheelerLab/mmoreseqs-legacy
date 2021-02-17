#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore-prep.sh	
#	AUTHOR:		David Rich
#	DESC: 		Creates any necessary files to prepare to run main pipeline.
###########################################################################

# PROGRAM:
# (1a) if <target> input is .msa file, convert to .hmm file
# (1b) if <target> input is .fasta file, convert to .hmm file
# (1c) verify that final <target> is .hmm file
#
# (2a) if <query> input is .hmm file, convert to .fasta file
# (2b) verify that final <query> is .fasta file
#
# (3a) if <target_mmseqs> is .msa file, convert to .mm_msa file.
# (3b) if <target_mmseqs> is .mm_msa file, convert to .hhm file.
# (3c) verify that final <target_mmseqs> file is .hhm or .fasta file.

# only echoes if has higher verbosity level
function echo_v
{
	if (( $VERBOSE >= $1 ))
	then
		echo $2
	fi
}
VERBOSE="${VERBOSE:-3}"

# finds the bench directory
function GET_BENCH_DIR
{
	local DIR="$(pwd)"
	echo ${DIR}
}
# get the script directory
BENCH_DIR=$(GET_BENCH_DIR)

# finds the script directory
function GET_SCRIPT_DIR
{
	local DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
	echo ${DIR}
}
SCRIPT_DIR_SIMPLE=$(GET_SCRIPT_DIR)

# commandline variables (verify proper number of variables)
NUM_ARGS=$#
if (( $NUM_ARGS < 3 )); then 
	echo "ERROR: illegal number of parameters: $NUM_ARGS of 3"
	echo "Usage: <tmp_root>"
	exit
fi

TMP_ROOT=$1

echo_v 1 "MMORE SCRIPT: $ARG_TARGET $ARG_QUERY $ARG_TARGET_MMSEQS"

BENCH_DIR=$(GET_BENCH_DIR)
SCRIPT_DIR=$(GET_SCRIPT_DIR_FULL)

source "${BENCH_DIR}/mmore-get-tools.sh"

# if args not supplied, falls back on these defaults
# main args:
TARGET="${TARGET:-$TARGET_ARG}"
QUERY="${QUERY:-$QUERY_ARG}"
TARGET_MMSEQS="${TARGET_MMSEQS:-$TARGET_MMSEQS_ARG}"
# main arg types:
TARGET_TYPE="${TARGET_TYPE:-'HMM'}"
QUERY_TYPE="${QUERY_TYPE:-'FASTA'}"
TARGET_MMSEQS_TYPE="${TARGET_MMSEQS_TYPE:-'HHM'}"
# options :
ROOT_DIR="${ROOT_DIR:-'/'}"
TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-XXXX)}"
RM_TEMP="${RM_TEMP:-0}" 
# interrim files:
# IS_HMM_OUT="${IS_HMM_OUT:-FALSE}"
# IS_MM_MSA_OUT="${IS_MM_MSA_OUT:-FALSE}"
# IS_HHM_OUT="${IS_HHM_OUT:-FALSE}"
# output files:

# ##### MAIN ###### #

# temporary directory
TMP_ROOT=${TEMP_DIR}/
TMP=${TMP_ROOT}/

# temporary subdirectories for mmseqs and 
TMP_MMSEQS=${TMP}/mmseqs/
TMP_=${TMP}/cloud/

# temporary databases for mmseqs
QUERY_MMSEQS=${TMP_MMSEQS}/query
TARGET_MMSEQS=${TMP_MMSEQS}/target

# result destination files
RESULT_RAW_MMSEQS=${TMP_MMSEQS}/results
RESULT_MMSEQS=${TMP_MMSEQS}/results.m8
RESULT_PLUS_MMSEQS=${TMP_MMSEQS}/results.m8+
RESULT_CLOUD=${TMP_CLOUD}/results.m8

# make tmp directories
mkdir $TMP_ROOT
mkdir $TMP
mkdir $TMP_CLOUD
mkdir $TMP_MMSEQS

# mmseqs options
MMSEQS_SHUFFLE=""

# === CLEAN UP === #
# remove temporary folders
if ( $REMOVE_TEMP > 0 )
then
	rm -r $TEMP_DIR
fi
