#!/usr/bin/bash
###########################################################################
#	 - FILE:  mmore_functions.sh	
#	 - DESC:  MMORE helper functions.
###########################################################################

# only echoes if has higher verbosity level
function echo_v # (verbose_level, message)
{
	local V="$1"
	local MESSAGE="$2"
	if (( $VERBOSE >= $V ))
	then
		echo "$MESSAGE"
	fi
}
VERBOSE="${VERBOSE:-3}"
if (( $VERBOSE > 3 )); then 
	VERB_ARG="-v"
fi

# finds the bench directory
function GET_BENCH_DIR
{
	local DIR="$(pwd)"
	echo ${DIR}
}
# get the script directory
BENCH_DIR=$(GET_BENCH_DIR)

# echo the command before running it
function ECHO_AND_RUN # (verbose_level, command_args...)
{
	# local V=$1

	# if command output location set, then append to it.
	if [[ -n $COMMAND_OUTPUT ]]
	then
		echo "${@:2}" >> $COMMAND_OUTPUT
	fi 

	# if verbosity is high enough, output command
	if (( $VERBOSE >= $1 ));
	then
		echo "# COMMAND: ${@:2}"
	fi

	# if performing dry run, dont run command
	if (( "PIPELINE_DRY_RUN" == "1" ))
	then 
		:
	else
		time "${@:2}"
	fi
}

# toggle 1<->
function TOGGLE # (bool)
{
	local BOOL=$1

	if (( $BOOL == 0 ));
	then
		echo "1"
		return
	elif (( $BOOL >= 1 ));
	then 
		echo "0"
		return
	fi

	echo "ERROR: Not boolean not correct value."
}

# finds the script directory
function GET_SCRIPT_DIR
{
	local DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
	echo ${DIR}
}
SCRIPT_DIR_SIMPLE=$(GET_SCRIPT_DIR)

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

# get time in seconds (with millisecond precision)
function GET_TIME
{
	echo $( date +%s.%N )
}

# get duration in seconds (with millisecond precision)
function GET_DURATION # (begin_time, end_time)
{
	local BEG_TIME=$1
	local END_TIME=$2
	local DUR_TIME="0"
	# DUR_TIME=$( echo "${END_TIME} - ${BEGIN_TIME}" | bc -l )
	DUR_TIME=$( python <<< "duration = ${END_TIME} - ${BEG_TIME}; print(duration)" ) 
	echo $DUR_TIME
}

# make directory safely
function MAKE_DIR # (folder_name)
{
	local DIR=$1
	if [[ ! -e $DIR ]]; then
		echo_v 3 "Making directory: '$DIR'..."
		mkdir $DIR
	elif [[ ! -d $DIR ]]; then
		echo "ERROR: '$DIR' already exists but is not a directory." 1>&2
	else
		:
		# echo_v 2 "Warning: Directory '$DIR' already exists."
	fi
}

# get file from full path
function GET_FILE_FROM_PATH
{
	local FULLFILE=$1
	local FILE=$(basename -- "$FULLFILE")
	echo $FILE
}

# get file extension
function GET_FILE_EXTENSION
{
	local FULLFILE=$1
	local FILE=$(GET_FILE_FROM_PATH "$FULLFILE")
	local EXT="${FILE##*.}"
	echo $EXT
}

# get base file name (w/o extension)
function GET_FILE_NAME
{
	local FULLFILE=$1
	local FILE=$(GET_FILE_FROM_PATH "$FULLFILE")
	local NAME="${FILE%*.}"
	echo $NAME
}

# infer file type based on file extension
function INFER_FILETYPE # (full_filename)
{
	local FULLFILE=$1
	local TYPE=""
	local EXT=$(GET_FILE_EXTENSION $FULLFILE)
	EXT=${EXT^^}
	if [[ $EXT == "HMM" ]]; then 
		TYPE="HMM"
	elif [[ $EXT == "FASTA" || $EXT == "FA" ]]; then 
		TYPE="FASTA"
	elif [[ $EXT == "HHM" ]]; then 
		TYPE="HHM"
	elif [[ $EXT == "MSA" ]]; then 
		TYPE="MSA"
	elif [[ $EXT == "MM_MSA" ]]; then 
		TYPE="MM_MSA"
	elif [[ $EXT == "MMDB" ]]; then 
		TYPE="MMDB"
	elif [[ $EXT == "MMDB_S" ]]; then 
		TYPE="MMDB_S"
	elif [[ $EXT == "MMDB_P" ]]; then 
		TYPE="MMDB_P"
	fi

	echo $TYPE
}

# check if file exists ()
function CHECK_FILE_EXISTS # (file)
{
	local FILE=$1
	local IF_FILE=""
	if [[ -f $FILE ]]; then
		IF_FILE="EXISTS"
	fi 
	echo $IF_FILE
}

# append argument with flag if not null
function SET_ARGOPT # (num_required_args, option_flag, option_args...)
{
	local REQARGS="${1}"
	local OPT="${2}"
	local ARG="${@:3}"
	local NUMARGS=$((${#} - 2))

	# if ARG is not null, append it 
	if (( $NUMARGS == $REQARGS ))
	then 
		echo "$OPT $ARG"
	else 
		echo "ERROR: Incorrect number of args (${NUMARGS} of ${REQARGS}: ${@})."
	fi
}

# if append argument with flag if not null
function IF_SET_ARGOPT  # (condition, num_required_args, option_flag, options_args...)
{
	local IF_SET="${1}"
	local ARGOPTS="${@:2}"

	# if ARG is not null, append it 
	if (( $IF_SET == 1 ))
	then
		SET_ARGOPT $ARGOPTS
	fi
}

# copies file, either through copying or symbolic link
function COPY_FILE # (source, destination)
{
	if [ "$DO_COPY" == "0" ] 
	then
		COPY="ln -s"
	else 
		COPY="cp"
	fi 

	local SRC="$1"
	local DEST="$2"

	$COPY "$SRC" "$DEST"
}

# if exit code is not zero (success), output exit code and terminate program
function CHECK_ERROR_CODE
{
	local ERROR_CODE="${EXIT_CODE:-"$?"}"
	local FUNCTION="${FUNCTION:-"$1"}"
	echo_v 3 "# ERROR_CODE: $ERROR_CODE"

	if (( EXIT_CODE != 0 ))
	then
		echo "ERROR: Executed ${FUNCTION} in ${PROGRAM} in unsuccessfully with ERROR_CODE:(${EXIT_CODE})."
		if (( DO_IGNORE_ERRORS == 0 )) 
		then 
			exit $EXIT_CODE
		fi
	fi 
}

# copy if condition is met 
function IF_THEN_COPY # (source, destination, conditional)
{
	local REQ_ARGS=3
	if (( $# < $REQ_ARGS ))
	then
		echo "WARNING: Did not copy. Insufficient Args: ($# of $REQ_ARGS)"
		return
	fi

	local SRC="$1"
	local DEST="$2"
	local COND="$3"

	if [[ -n $COND ]] && (( $COND != 0 ))
	then
		cp "$SRC" "$DEST"
	else
		echo_v 3 "WARNING: Did not copy ($SRC -> $DEST). Condition not met (COND = $COND)."
	fi
}

# convert eval to pval
function PVAL_TO_EVAL # (pval, dbsize)
{
	local PVAL=$1
	local DBSIZE=$2

	if (( $# != 2 )); then 
		echo "#ERROR: Too few arguments."
	fi 

	if [[ $DBSIZE == "inf" ]]; then
		DBSIZE="np.inf"
	elif [[ $DBSIZE == "-inf" ]]; then 
		DBSIZE="-np.inf"
	fi

	local EVAL=$( python <<< "import numpy as np; eval = ${PVAL} * ${DBSIZE}; print('%.3e' % eval)" )

	echo $EVAL
}

# set default parameters for uninitialized environmental variables
function SET_ENV_ARG_DEFAULTS
{
	# process environmental variables: if optional args not supplied by environment, falls back on these defaults
	{
		# pipeline:
		# pipeline options:
		PIPELINE_DRY_RUN="${PIPELINE_DRY_RUN:-"0"}"
		COMMAND_OUTPUT="${COMMAND_OUTPUT:-""}"
		# tool paths:
		MMORESEQS="${MMORESEQS:-"mmoreseqs"}"
		MMSEQS="${MMSEQS:-"mmseqs"}"
		HMMSEARCH="${HMMSEARCH:-"hmmsearch"}"
		HMMBUILD="${HMMBUILD:-"hmmbuild"}"
		HMMEMIT="${HMMEMIT:-"hmmemit"}"
		# mmore args:
		TARGET_MMORE_P="${TARGET_MMORE_P:-""}"
		TARGET_MMORE_S="${TARGET_MMORE_S:-""}"
		QUERY_MMORE="${QUERY_MMORE:-""}"
		TARGET_MMORE_P_TYPE="${TARGET_MMORE_P_TYPE:-""}"
		TARGET_MMORE_S_TYPE="${TARGET_MMORE_S_TYPE:-""}"
		QUERY_MMORE_TYPE="${QUERY_MMORE_TYPE:-""}"
		# mmseqs args:
		TARGET_MMSEQS_P="${TARGET_MMSEQS_P:-""}"
		TARGET_MMSEQS_S="${TARGET_MMSEQS_S:-""}"
		QUERY_MMSEQS="${QUERY_MMSEQS:-""}"
		TARGET_MMSEQS_P_TYPE="${TARGET_MMSEQS_P_TYPE:-""}"
		TARGET_MMSEQS_S_TYPE="${TARGET_MMSEQS_S_TYPE:-""}"
		QUERY_MMSEQS_TYPE="${QUERY_MMSEQS_TYPE:-""}"
		# prep args:
		PREP_DIR="${PREP_DIR:-""}"
		TARGET_PREP="${TARGET_PREP:-""}"
		QUERY_PREP="${QUERY_PREP:-""}"
		TARGET_PREP_TYPE="${TARGET_PREP_TYPE:-"MSA"}"
		QUERY_PREP_TYPE="${QUERY_PREP_TYPE:-"FASTA"}"
		# other input args:
		TARGET_FASTA="${TARGET_FASTA:-""}"
		TARGET_HMM="${TARGET_HMM:-""}"
		QUERY_FASTA="${QUERY_FASTA:-""}"
		# data input:
		TARGET_DBSIZE="${TARGET_DBSIZE:-""}"
		QUERY_DBSIZE="${QUERY_DBSIZE:-""}"
		PREFILTER_DBSIZE="${PREFILTER_DBSIZE:-""}"
		P2S_DBSIZE="${P2S_DBSIZE:-""}"
		S2S_DBSIZE="${S2S_DBSIZE:-""}"
		MMORE_DBSIZE="${MMORE_DBSIZE:-""}"
		# time input:
		TIME_PREP="${TIME_PREP:-"0"}"
		TIME_MMSEQS="${TIME_MMSEQS:-"0"}"
		TIME_MMSEQS_PREFILTER="${TIME_MMSEQS_PREFILTER:-"0"}"
		TIME_MMSEQS_P2S="${TIME_MMSEQS_P2S:-"0"}"
		TIME_MMSEQS_IDS="${TIME_MMSEQS_IDS:-"0"}"
		TIME_MMSEQS_S2S="${TIME_MMSEQS_S2S:-"0"}"
		TIME_MMSEQS_CONVERTALIS="${TIME_MMSEQS_CONVERTALIS:-"0"}"
		TIME_MMORE="${TIME_MMORE:-"0"}"
		# directories:
		ROOT_DIR="${ROOT_DIR:-"/"}"
		SCRIPT_DIR="${SCRIPT_DIR:-""}"
		TEMP_DIR="${TEMP_DIR:-""}"
		# pipeline options:
		RM_TEMP="${DO_RM_TEMP:-"0"}" 
		DO_PREP="${DO_PREP:-"0"}"
		DO_COPY="${DO_COPY:-"1"}"
		DO_STATS="${DO_STATS:-"0"}"
		DO_OVERWRITE="${DO_OVERWRITE:-"1"}"
		DO_IGNORE_WARNINGS="${DO_IGNORE_WARNINGS:-"1"}"
		DO_IGNORE_ERRORS="${DO_IGNORE_ERRORS:-"0"}"
		NUM_THREADS="${NUM_THREADS:-"1"}"
		SEARCH_TYPE="${SEARCH_TYPE:-"P2S"}"

		# mmseqs parameters:
		MMSEQS_MAXSEQS="${MMSEQS_MAXSEQS:-"1000"}"
		MMSEQS_PVAL="${MMSEQS_PVAL:-"1e-2"}"
		MMSEQS_EVAL="${MMSEQS_EVAL:-""}"
		# mmseqs options:
		MMSEQS_DO_MMSEQS="${MMSEQS_DO_MMSEQS:-"1"}"
		MMSEQS_NUM_THREADS="${MMSEQS_NUM_THREADS:-"1"}"
		MMSEQS_DO_VIT_ALN="${MMSEQS_DO_VIT_ALN:-"0"}"		
		# mmseqs prefilter parameters:
		MMSEQS_PREFILTER_MAXSEQS="${MMSEQS_PREFILTER_MAXSEQS:-"1000"}"
		MMSEQS_KMER="${MMSEQS_KMER:-"7"}"
		MMSEQS_KSCORE="${MMSEQS_KSCORE:-"80"}"
		MMSEQS_DO_UNGAPPED="${MMSEQS_DO_UNGAPPED:-"1"}"
		MMSEQS_UNGAPPED="${MMSEQS_UNGAPPED:-"15"}"
		MMSEQS_DO_SENS="${MMSEQS_DO_SENS:-"0"}"
		MMSEQS_SENS="${MMSEQS_SENS:-"7.5"}"
		MMSEQS_ALTALIS="${MMSEQS_ALTALIS:-"0"}"
		# mmseqs p2s parameters:
		MMSEQS_P2S_ALTALIS="${MMSEQS_P2S_ALTALIS:-"0"}"
		MMSEQS_P2S_PVAL="${MMSEQS_P2S_PVAL:-"1e-2"}"
		MMSEQS_P2S_EVAL="${MMSEQS_P2S_PVAL:-""}"
		# mmseqs s2s parameters:
		MMSEQS_S2S_ALTALIS="${MMSEQS_S2S_ALTALIS:-"0"}"
		MMSEQS_S2S_PVAL="${MMSEQS_S2S_PVAL:-"inf"}"
		MMSEQS_S2S_EVAL="${MMSEQS_S2S_PVAL:-""}"
		# convert mmseqs to hmmer models
		MMSEQS_DO_CONVERT="${MMSEQS_DO_CONVERT:-"1"}"

		# mmore parameters:
		MMORE_ALPHA="${MMORE_ALPHA:-"12.0"}"
		MMORE_BETA="${MMORE_BETA:-"16.0"}"
		MMORE_GAMMA="${MMORE_GAMMA:-"5"}"
		MMORE_VITERBI_PVAL="${MMORE_VITERBI_PVAL:-"1e-3"}"
		MMORE_CLOUD_PVAL="${MMORE_CLOUD_PVAL:-"1e-5"}"
		MMORE_BOUNDFWD_PVAL="${MMORE_BOUNDFWD_PVAL:-"1e-5"}"
		MMORE_VITERBI_EVAL="${MMORE_VITERBI_EVAL:-"1e-3"}"
		MMORE_CLOUD_EVAL="${MMORE_CLOUD_EVAL:-"1e-5"}"
		MMORE_BOUNDFWD_EVAL="${MMORE_BOUNDFWD_EVAL:-"1e-5"}"
		# mmore options:
		MMORE_DO_MMORE="${MMORE_DO_MMORE:-"1"}"
		MMORE_NUM_THREADS="${MMORE_NUM_THREADS:-"1"}"
		MMORE_DO_FILTER="${MMORE_DO_FILTER:-"0"}"
		MMORE_DO_VIT_FILTER="${MMORE_DO_VIT_FILTER:-"0"}"
		MMORE_DO_CLD_FILTER="${MMORE_DO_CLD_FILTER:-"0"}"
		MMORE_DO_FWD_FILTER="${MMORE_DO_FWD_FILTER:-"0"}"
		MMORE_DO_BIAS="${MMORE_DO_BIAS:-"1"}"
		MMORE_DO_DOMAIN="${MMORE_DO_DOMAIN:-"1"}"
		MMORE_DO_FULL="${MMORE_DO_FULL:-"0"}"
		MMORE_DO_ALLOUT="${MMORE_DO_ALLOUT:-"1"}"
		MMORE_DO_MMSEQS_ALN="${MMORE_DO_MMSEQS_ALN:-"0"}"
		MMORE_DO_VIT_ALN="${MMORE_DO_VIT_ALN:-"0"}"
		MMORE_DO_POST_ALN="${MMORE_DO_POST_ALN:-"1"}"
		# interrim files (input):
		TARGET_MMSEQS_S_IN="${TARGET_MMSEQS_S_IN:-""}"
		TARGET_MMSEQS_P_IN="${TARGET_MMSEQS_P_IN:-""}"
		QUERY_MMSEQS_IN="${QUERY_MMSEQS_IN:-""}"
		TARGET_INDEX_IN="${TARGET_INDEX_IN:-""}"
		QUERY_INDEX_IN="${QUERY_INDEX_IN:-""}"
		MMSEQS_P2S_IN="${MMSEQS_P2S_IN:-""}"
		MMSEQS_S2S_IN="${MMSEQS_S2S_IN:-""}"
		MMSEQS_M8_IN="${MMSEQS_M8_IN:-""}"
		# interrim files (output):
		TARGET_MMSEQS_S_OUT="${TARGET_MMSEQS_S_OUT:-""}"
		TARGET_MMSEQS_P_OUT="${TARGET_MMSEQS_P_OUT:-""}"
		QUERY_MMSEQS_OUT="${QUERY_MMSEQS_OUT:-""}"
		TARGET_INDEX_OUT="${TARGET_INDEX_OUT:-""}"
		QUERY_INDEX_OUT="${QUERY_INDEX_OUT:-""}"
		MMSEQS_P2S_OUT="${MMSEQS_P2S_OUT:-""}"
		MMSEQS_S2S_OUT="${MMSEQS_S2S_OUT:-""}"
		MMSEQS_M8OUT_OUT="${MMSEQS_M8_OUT:-""}"
		# output files:
		MMORE_STDOUT_OUT="${MMORE_STDOUT_OUT:-""}"
		MMORE_STDERR_OUT="${MMORE_STDERR_OUT:-""}"
		MMORE_M8OUT_OUT="${MMORE_M8OUT_OUT:-""}"
		MMORE_MYOUT_OUT="${MMORE_MYOUT_OUT:-""}"
		MMORE_MYDOM_OUT="${MMORE_MYDOM_OUT:-""}"
		MMORE_MYTIME_OUT="${MMORE_MYTIME_OUT:-""}"
		MMORE_MYTHRESH_OUT="${MMORE_MYTHRESH_OUT:-""}"
	}
}

# set default parameters for uninitialized environmental variables
function SET_TOOLS
{
	echo "# SET TOOLS:"
	echo "  MMORESEQS || $MMORESEQS" 
	echo "            || $(which $MMORESEQS)"
	echo "   HMMBUILD || $HMMBUILD" 
	echo "            || $(which $HMMBUILD)"
	echo "     MMSEQS || $MMSEQS" 
	echo "            || $(which $MMSEQS)"
}

echo_v 3 "# mmore-functions.sh Loaded."
