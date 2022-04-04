#!/usr/bin/bash
###########################################################################
#	 - FILE:  mmore-prepsearch.sh	
#	 - DESC:   Run search with prepped file as input.
###########################################################################

# PROGRAM:
# (1) Extract Files from Prep Folder
# (2) Run Search

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
		echo "Error Message: $ERROR."

		exit $EXIT_CODE
	}
}

# ##### MAIN ###### #
{
	echo_v 2 "Running Program: mmore-prepsearch.sh..."
	PROGRAM+="> mmore-prepsearch.sh >"

	# load script dependencies 
	{
		# load functions
		echo_v 3 "Importing 'mmore_functions.sh..."
		LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_functions.sh"
	}

	BEG_TIME_PREPSEARCH=$(GET_TIME)

	# variable preprocessing 
	{
		# process commandline variables
		{
			#verify proper number of variables
			NUM_ARGS=$#
			REQ_ARGS=1
			if (( $NUM_ARGS < $REQ_ARGS )); then 
				echo "ERROR: Illegal number of main args: ($NUM_ARGS of $REQ_ARGS)"
				echo "Usage: <i:prep_folder>"
				echo "Accepted Filetypes: FASTA, MSA"
				exit 1
			fi

			# main commandline args
			ARG_PREP_DIR=$1
		}

		# process environmental variables: if optional args not supplied by environment, falls back on these defaults
		{
			# main args:
			PREP_DIR="${PREP_DIR:-$ARG_PREP_DIR}"
			TEMP_DIR="${PREP_DIR:-$ARG_PREP_DIR}"

			# extract files from prep directory
			TARGET_MMORE="${PREP_DIR}/mmore/db/target.hmm"
			TARGET_MMORE_P="${PREP_DIR}/mmore/db/target.hmm"
			TARGET_MMORE_S="${PREP_DIR}/mmore/db/target.fasta"
			QUERY_MMORE="${PREP_DIR}/mmore/db/query.fasta"
			TARGET_MMSEQS_P="${PREP_DIR}/mmseqs/db/target.pmmdb"
			TARGET_MMSEQS_S="${PREP_DIR}/mmseqs/db/target.smmdb"
			QUERY_MMSEQS="${PREP_DIR}/mmseqs/db/query.smmdb"
			COMMAND_OUTPUT="${PREP_DIR}/commands.txt"

			# default args:
			SET_ENV_ARG_DEFAULTS
			SET_TOOLS
		}

		# verify files exist 
		{
			:
		}

		# report variables
		{
			if (( $VERBOSE >= 3 )); then
				echo "# ============ MMORESEQS: PREPSEARCH ============"
				echo "# PREP_FOLDER:  $PREP_DIR"
			fi
		}
	}

	# Run MMORESEQS Search
	{
		# Call C Program
		# $MMORE mmore-search 									\

		# Call Bash Script
		LOAD_SOURCE "${SCRIPT_DIR}/mmoreseqs-search.sh" 	\
		"$TARGET_MMORE" 												\
		"$QUERY_MMORE" 												\
		"$TARGET_MMSEQS_P"  "$TARGET_MMSEQS_S" 				\
		"$QUERY_MMSEQS" 												\
		"$PREP_DIR" 													\

		CHECK_ERROR_CODE "MMORE_SEARCH"
	}

	END_TIME_PREPSEARCH=$(GET_TIME)
	TIME_PREP=$(GET_DURATION $BEG_TIME_PREPSEARCH $END_TIME_PREPSEARCH)
	printf "# TIME_PREP_SEARCH: %.3f sec\n" $TIME_PREP
}
