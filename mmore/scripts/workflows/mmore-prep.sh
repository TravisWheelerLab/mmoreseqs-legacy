#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore-prep.sh	
#	AUTHOR:		David Rich
#	DESC: 		Creates any necessary files to prepare to run main pipeline.
###########################################################################

# PROGRAM:
# (1a) if <target> is .msa file, convert to .hmm file
# (1b) if <target> is .fasta file, convert to .hmm file
# (1c) verify that final <target> is .hmm file
#
# (2a) if <query> is .msa file, convert to .hmm file
# (2b) if <query> is .hmm file, convert to .fasta (consensus) file
# (2c) verify that final <query> is .fasta file
#
# (3a) if <target_mmseqs> is .msa file, convert to .mm_msa file.
# (3b) if <target_mmseqs> is .mm_msa file, convert to .hhm file.
# (3c) if <target_mmseqs> is .hhm file, convert to .mm_db file.
# (3d) if <target_mmseqs> is .fasta file, convert to .mm_db file.
# (3e) verify that final <target_mmseqs> file is .mm_db file.
# 
# (4a) if <query_mmseqs> is .msa file, convert to .fasta (consensus) file.
# (4b) if <query_mmseqs> is .mm_msa file, convert to .fasta (consensus) file (?).
# (4c) if <query_mmseqs> is .hhm file, convert to .fasta (consensus) file (?).
# (4d) if <query_mmseqs> is .fasta file, convert to .mm_db file.
# (4e) verify that final <query_mmseqs> file is .mm_db file.
#
# (5a) 
# (5b)

# ##### FUNCTIONS ###### #
# main script only contains the minimum functions needed to load helper function script
{
	#set verbosity
	VERBOSE="${VERBOSE:-3}"
	if (( $VERBOSE > 3 )); then 
		VERB_ARG="-v"
	fi

	# only echoes if has higher verbosity level
	function echo_v
	{
		if (( $VERBOSE >= $1 ))
		then
			echo $2
		fi
	}

	# check if file exists ()
	function CHECK_FILE_EXISTS
	{
		local FILE=$1
		local FILE_EXISTS=""
		if [[ -f $FILE ]]; then
			FILE_EXISTS="TRUE"
		fi 
		echo $FILE_EXISTS
	}

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
	SCRIPT_DIR=$(GET_SCRIPT_DIR)

	# finds the script directory (across source, aliasing, and symlinks)
	function GET_SCRIPT_DIR_FULL
	{
		local SOURCE="${BASH_SOURCE[0]}"
		while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
			local DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
			SOURCE="$(readlink "$SOURCE")"
			[[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
		done
		DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
		echo ${DIR}
	}
	# SCRIPT_DIR=$(GET_SCIPT_DIR_FULL)

	# verifies that source file exists; if so, then imports source
	function LOAD_SOURCE 
	{
		local SOURCE=$1
		local FILE_EXISTS=$(CHECK_FILE_EXISTS $SOURCE)
		if [ -z FILE_EXISTS ]; then
			echo "ERROR: Source file '$SOURCE' does not exist."
		fi 
		source $SOURCE "${@:2}"
	}

	# in the event of an error, when want to clean up working files
	function ERROR_EXIT 
	{
		local EXIT_CODE=$1
		echo "Exiting program with error code: $EXIT_CODE."

		exit $EXIT_CODE
	}
}

# ##### MAIN ###### #
echo_v 1 "Running Program: mmore-prep.sh..."

# load functions
echo_v 3 "Importing 'mmore_functions.sh..."
LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_functions.sh"
# load tools 
echo_v 3 "Importing 'mmore_get-tools.sh'..."
LOAD_SOURCE "${SCRIPT_DIR}/mmore_get-tools.sh"

# variable preprocessing 
{
	# process commandline variables
	{
		#verify proper number of variables
		NUM_ARGS=$#
		REQ_ARGS=4
		if (( $NUM_ARGS < $REQ_ARGS )); then 
			echo "ERROR: Illegal number of main args: ($NUM_ARGS of $REQ_ARGS)"
			echo "Usage: <target> <query> <target_mmseqs> <query_mmseqs> | <target_type> <query_type> <target_mmseqs_type> <query_mmseqs_type>"
			echo "Accepted Filetypes: FASTA, HMM, HHM, MSA, MM_MSA, MM_DB"
			exit 1
		fi

		# main commandline args
		ARG_TARGET=$1
		ARG_QUERY=$2
		ARG_TARGET_MMSEQS=$3
		ARG_QUERY_MMSEQS=$4
		# optional commandline args
		ARG_TARGET_TYPE=$5
		ARG_QUERY_TYPE=$6
		ARG_TARGET_MMSEQS_TYPE=$7
		ARG_QUERY_MMSEQS_TYPE=$8

		# print command line arguments
		echo_v 1 "mmore.sh: [0]$ARG_TARGET [1]$ARG_QUERY [2]$ARG_TARGET_MMSEQS [3]$ARG_QUERY_MMSEQS | [4]$ARG_TARGET_TYPE [5]$ARG_QUERY_TYPE [6]$ARG_TARGET_MMSEQS_TYPE [7]$ARG_QUERY_MMSEQS_TYPE"
	}

	# process environmental variables: if optional args not supplied by environment, falls back on these defaults
	{
		# main args:
		TARGET="${TARGET:-$ARG_TARGET}"
		QUERY="${QUERY:-$ARG_QUERY}"
		TARGET_MMSEQS="${TARGET_MMSEQS:-$ARG_TARGET_MMSEQS}"
		QUERY_MMSEQS="${TARGET_MMSEQS:-$ARG_TARGET_MMSEQS}"
		# main arg types:
		TARGET_TYPE="${TARGET_TYPE:-$ARG_TARGET_TYPE}"
		QUERY_TYPE="${QUERY_TYPE:-$ARG_QUERY_TYPE}"
		TARGET_MMSEQS_TYPE="${TARGET_MMSEQS_TYPE:-$ARG_TARGET_MMSEQS_TYPE}"
		QUERY_MMSEQS_TYPE="${QUERY_MMSEQS_TYPE:-$ARG_QUERY_MMSEQS_TYPE}"
		# options :
		ROOT_DIR="${ROOT_DIR:-'/'}"
		TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-XXXX)}"
		RM_TEMP="${RM_TEMP:-0}" 
		DO_PREP="${DO_PREP:-}"
		DO_OVERWRITE="${DO_OVERWRITE:-1}"
		DO_IGNORE_WARNINGS="${DO_IGNORE_WARNINGS:-1}"
		NUM_THREADS="${NUM_THREADS:-1}"
		# mmseqs parameters:
		KMER="${KMER:-7}"
		K_SCORE="${K_SCORE:-80}"
		MIN_UNGAPPED_SCORE="${UNGAPPEDVIT_THRESH:-15}"
		PVAL_CUTOFF="${GAPPEDVIT_THRESH:-0.001}"
		EVAL_CUTOFF="${EVAL_CUTOFF:-200}"
		SENSITIVITY="${SENSITIVITY:-7.5}"
		# mmore parameters:
		ALPHA="${ALPHA:-12.0}"
		BETA="${BETA:-16.0}"
		GAMMA="${GAMMA:-5}"
		VIT_THRESH="${VIT_THRESH:-1e-3}"
		CLOUD_THRESH="${CLOUD_THRESH:-1e-5}"
		FWD_THRESH="${FWD_THRESH:-1e-5}"
		DO_FILTER="${DO_FILTER:-0}"
		DO_BIAS="${DO_BIAS:-1}"
		DO_DOMAIN="${DO_DOMAIN:-1}"
		DO_FULL="${DO_FULL:-0}"
		# interrim files:
		TARGET_MMSEQS_DB="${TARGET_MMSEQS_DB:-""}"
		QUERY_MMSEQS_DB="${QUERY_MMSEQS_DB:-""}"
		TARGET_INDEX="${TARGET_INDEX:-""}"
		QUERY_INDEX="${QUERY_INDEX:-""}"
		MMSEQS_P2SOUT="${MMSEQS_P2SOUT:-""}"
		MMSEQS_S2SOUT="${MMSEQS_S2SOUT:-""}"
		MMSEQS_M8OUT="${MMSEQS_OUT:-""}"
		# output files:
		MMORE_STDOUT="${MMORE_STDOUT:-""}"
		MMORE_M8OUT="${MMORE_M8OUT:-""}"
		MMORE_MYOUT="${MMORE_MYOUT:-""}"
		MMORE_TIMEOUT="${MMORE_TIMEOUT:-""}"
	}

	# list of all main arguments
	{
		declare -A MAIN_FILES
		MAIN_FILES["TARGET"]="$TARGET"
		MAIN_FILES["QUERY"]="$QUERY" 
		MAIN_FILES["TARGET_MMSEQS"]="$TARGET_MMSEQS" 
		MAIN_FILES["QUERY_MMSEQS"]="$QUERY_MMSEQS"

		declare -A MAIN_TYPES
		MAIN_TYPES["TARGET_TYPE"]="$TARGET_TYPE"
		MAIN_TYPES["QUERY_TYPE"]="$QUERY_TYPE" 
		MAIN_TYPES["TARGET_MMSEQS_TYPE"]="$TARGET_MMSEQS_TYPE" 
		MAIN_TYPES["QUERY_MMSEQS_TYPE"]="$QUERY_MMSEQS_TYPE"
	}

	# list of all optional arguments
	{
		declare -A ARG_OPT
	}

	# specify input files (will verify they exist during preprocessing)
	{
		declare -A INPUT_FILES
	}

	# specify temporary files (will delete during cleanup stage)
	{
		declare -A TEMP_FILES
	}

	# specify output files (will move/copy from temp folder during cleanup stage)
	{
		OUTPUT_FILES=""
	}

	# report variables
	{
		if (( $VERBOSE >= 3 )); then
			echo "#        TARGET: $TARGET"
			echo "#         QUERY: $QUERY"
			echo "#       RESULTS: $RESULTS"
			
			echo "#  FILTER_SCORE: $K_SCORE"
			echo "#   UNGAP_SCORE: $MIN_UNGAPPED_SCORE"
			echo "#     GAP_SCORE: $PVAL_CUTOFF => $EVAL_CUTOFF"
			
			echo "#         ALPHA: $ALPHA"
			echo "#          BETA: $BETA"
			echo "#         GAMMA: $GAMMA"
			echo "#  REPORT_SCORE: $PVALUE_REPORT => $EVALUE_REPORT"
		fi
	}

	# build temporary folders
	{
		# top-level temporary directory
		TMP_ROOT=${TEMP_DIR}/
		TMP=${TMP_ROOT}/
		# temporary subdirectories for mmseqs and mmore
		TMP_MMSEQS=${TMP}/mmseqs/
		TMP_MMORE=${TMP}/mmore/
		# temporary subdirectories for mmseqs databases
		TMP_MMSEQS_DB=${TMP_MMSEQS}/db/

		# make tmp directories
		MAKE_DIR $TMP_ROOT
		MAKE_DIR $TMP_MMORE
		MAKE_DIR $TMP_MMSEQS
		MAKE_DIR $TMP_MMSEQS_DB

		# list of temp folders 
		TEMP_FOLDERS=""
		TEMP_FOLDERS+="${TMP}" 
		TEMP_FOLDERS+="${TMP_MMSEQS}" 
		TEMP_FOLDERS+="${TMP_MMORE}"
		TEMP_FOLDERS+="${TMP_MMSEQS_DB}"
	}
}

# === FILE FORMATTING === #
{
	echo "# (0/4) File formatting and building..."

	# QUERY file formatting
	{
		# currently not supported
		if [ "$QUERY_TYPE" == "MSA" ]
		then
			QUERY_MSA=$QUERY
			echo "WARNING: QUERY is an MSA file. Conversion currently not supported."	
			exit 1	
		fi

		# If query is a hmm file, then convert to consensus fasta file
		# WARNING: this may not be intended
		if [[ "$QUERY_TYPE" = "HMM" ]]
		then
			QUERY_HMM=$QUERY
			echo_v 1 	"WARNING: QUERY is an HMM file. Profile-to-Profile is not supported, \
							so query being converted to FASTA consensus sequences. 				  \
							This will result in a loss of accuracy."			

			QUERY_TYPE="FASTA"
		fi

		# if target mmseqs file is a hhm, then we are ready.
		if [[ "$QUERY_TYPE" = "FASTA" ]]
		then
			QUERY_FASTA=$QUERY
			echo_v 1 ""
			QUERY_TYPE="FASTA"
		fi

		# Final check if file is a valid file format
		if [ "$QUERY_TYPE" != "FASTA" ]
		then
			echo "ERROR: QUERY is an unsupported filetype of '$QUERY_TYPE'."
			echo "Valid QUERY Types: FASTA, HMM, MSA."
			exit 1
		fi
	}

	# TARGET file formatting
	{
		# if target is a msa file, then convert to HMM
		if [[ "TARGET_TYPE" = "MSA" ]]
		then 
			TARGET_MSA=$TARGET 
			echo "WARNING: TARGET is an MSA file. Conversion currently not supported."
			exit 1
			TARGET_TYPE="MSA"
		fi

		# If target is a hmm file, then 
		if [[ "$TARGET_TYPE" != "HMM" ]]
		then
			TARGET_HMM=$TARGET
		fi

		# If target is a fasta file, then we need to 
		# build single sequence HMMs for all sequences in database. 
		if [[ "$TARGET_TYPE" = "FASTA" ]]
		then
			TARGET_FASTA=$TARGET 
		fi

		# Final check if file is a valid file format.
		if [[ "$TARGET_TYPE" != "HMM" ]]
		then
			echo "ERROR: Target is an unsupported type of '$TARGET_TYPE'."
			exit 1
		fi
	}

	# TARGET_MMSEQS file formatting
	{
		# if target mmseqs file is a msa file, then we need to
		# convert it to mm_msa, then generate 
		if [[ "$TARGET_MMSEQS_TYPE" == "MSA" ]]
		then
			TARGET_MMSEQS_MSA=$TARGET_MMSEQS

			$MMSEQS 

			TARGET_MMSEQS=$TARGET_MMSEQS_MM_MSA
			TARGET_MMSEQS_TYPE="MM_MSA"
		fi

		# if target mmseqs file is a msa file, then we need to convert it to mm_msa
		if [[ "$TARGET_MMSEQS_TYPE" == "MM_MSA" ]]
			TARGET_MMSEQS_MM_MSA=$TARGET_MMSEQS


		fi

		# if target mmseqs file is a hhm, then we are ready.
		if [[ "$TARGET_MMSEQS_TYPE" == "HHM" ]]
			TARGET_MMSEQS_HHM=$TARGET_MMSEQS
		fi

		# Final check if file is a valid file format
		if [ "$TARGET_MMSEQS_TYPE" != "MM_DB" ]
			echo "ERROR: MMseqs target is an unsupported type of '$TARGET_MMSEQS_TYPE'."
			exit 1
		fi
	}



}

