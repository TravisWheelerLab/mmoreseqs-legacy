#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore.sh	
#	AUTHOR:		David Rich
#	DESC: 		MMORE helper functions.
###########################################################################

# only echoes if has higher verbosity level
function echo_v
{
	if (( $VERBOSE >= $1 ))
	then
		echo $2
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
	echo date +%s.%N
}

# get duration in seconds (with millisecond precision)
function GET_DURATION
{
	local BEG_TIME=$1
	local END_TIME=$2
	echo $( echo "${END_TIME} - ${BEGIN_TIME}" | bc -l )
}

# make directory safely
function MAKE_DIR
{
	local DIR=$1
	if [[ ! -e $DIR ]]; then
		echo_v 3 "Making directory: '$DIR'..."
		mkdir $DIR
	elif [[ ! -d $DIR ]]; then
		echo "ERROR: '$DIR' already exists but is not a directory." 1>&2
	else
		echo_v 2 "Warning: Directory '$DIR' already exists."
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
function INFER_FILETYPE
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
function CHECK_FILE_EXISTS
{
	local FILE=$1
	local IF_FILE=""
	if [[ -f $FILE ]]; then
		IF_FILE="EXISTS"
	fi 
	echo $IF_FILE
}

# append argument with flag if not null
function ADD_OPT
{
	ARGOPTS=$1
	local OPT=$1
	local ARG=$2

	# if ARG is not null, append it 
	if [ -z $ARG ]
	then 
	ARGOPTS+=" $OPT $ARG "
	fi
}

# copies file, either through copying or symbolic link
function COPY_FILE
{
	if [ "$DO_COPY" == "0" ] 
	then
		COPY="ln -s"
	else 
		COPY="cp"
	fi 

	local SRC=$1
	local DEST=$2

	$COPY "$SRC" "$DEST"
}

# set default parameters for uninitialized environmental variables
function SET_ENV_ARG_DEFAULTS
{
	# process environmental variables: if optional args not supplied by environment, falls back on these defaults
	{
		# mmore args:
		TARGET_MMORE="${TARGET_MMORE:-""}"
		QUERY_MMORE="${QUERY_MMORE:-""}"
		TARGET_MMORE_TYPE="${TARGET_MMORE_TYPE:-""}"
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
		TIME_PREP="${TIME_PREP:-""}"
		TIME_MMSEQS="${TIME_MMSEQS:-""}"
		TIME_MMSEQS_PREFILTER="${TIME_MMSEQS_PREFILTER:-""}"
		TIME_MMSEQS_P2S="${TIME_MMSEQS_P2S:-""}"
		TIME_MMSEQS_S2S="${TIME_MMSEQS_S2S:-""}"
		TIME_MMORE=""
		# directories:
		ROOT_DIR="${ROOT_DIR:-'/'}"
		SCRIPT_DIR="${SCRIPT_DIR:-""}"
		TEMP_DIR="${TEMP_DIR:-""}"
		# pipeline options:
		RM_TEMP="${DO_RM_TEMP:-0}" 
		DO_PREP="${DO_PREP:-0}"
		DO_COPY="${DO_COPY:-1}"
		DO_OVERWRITE="${DO_OVERWRITE:-1}"
		DO_IGNORE_WARNINGS="${DO_IGNORE_WARNINGS:-1}"
		DO_IGNORE_ERRORS="${DO_IGNORE_ERRORS:-0}"
		NUM_THREADS="${NUM_THREADS:-1}"
		SEARCH_TYPE="${SEARCH_TYPE:-"P2S"}"
		# mmseqs parameters:
		MMSEQS_NUM_THREADS="${MMSEQS_NUM_THREADS:-1}"
		MMSEQS_PREFILTER_MAXSEQS="${MMSEQS_PREFILTER_MAXSEQS:-300}"
		MMSEQS_KMER="${MMSEQS_KMER:-7}"
		MMSEQS_KSCORE="${MMSEQS_KSCORE:-80}"
		MMSEQS_UNGAPPED="${MMSEQS_UNGAPPED:-15}"
		# MMSEQS_SENS="${MMSEQS_SENS:-7.5}"
		MMSEQS_S2S_MAXSEQS="${MMSEQS_P2S_MAXSEQS:-0}"
		MMSEQS_S2S_MAXSEQS="${MMSEQS_S2S_MAXSEQS:-5}"
		MMSEQS_PVAL="${MMSEQS_PVAL:-1e-3}"
		MMSEQS_PVAL="${MMSEQS_PVAL:-1e-1}"
		MMSEQS_EVAL="${MMSEQS_EVAL:-1e-3}"
		MMSEQS_PVAL="${MMSEQS_PVAL:-1e-1}"
		# mmseqs options:
		MMSEQS_DO_VIT_ALN="${MMSEQS_DO_VIT_ALN:-0}"		

		# mmore parameters:
		MMSEQS_NUM_THREADS="${MMSEQS_NUM_THREADS:-1}"
		MMORE_ALPHA="${MMORE_ALPHA:-12.0}"
		MMORE_BETA="${MMORE_BETA:-16.0}"
		MMORE_GAMMA="${MMORE_GAMMA:-5}"
		MMORE_VITERBI="${MMORE_VITERBI:-1e-3}"
		MMORE_CLOUD="${MMORE_CLOUD:-1e-5}"
		MMORE_BOUNDFWD="${MMORE_BOUNDFWD:-1e-5}"
		# mmore options:
		MMORE_DO_FILTER="${MMORE_DO_FILTER:-0}"
		MMORE_DO_BIAS="${MMORE_DO_BIAS:-1}"
		MMORE_DO_DOMAIN="${MMORE_DO_DOMAIN:-1}"
		MMORE_DO_FULL="${MMORE_DO_FULL:-0}"
		MMORE_DO_ALLOUT="${MMORE_DO_ALLOUT:-1}"
		MMORE_DO_MMSEQS_ALN="${MMORE_DO_MMSEQS_ALN:-0}"
		MMORE_DO_VIT_ALN='${MMORE_DO_VIT_ALN:-0}'
		MMORE_DO_POST_ALN='${MMORE_DO_POST_ALN:-1}'
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
		MMSEQS_M8_OUT="${MMSEQS_M8_OUT:-""}"
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

echo_v 3 "# mmore-functions.sh Loaded."
