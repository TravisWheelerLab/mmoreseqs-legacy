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
echo "VERBOSE LEVEL: $VERBOSE"

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

# get time in milliseconds
function GET_TIME
{
	echo date +%s.%N
}

# make directory safely
function MAKE_DIR
{
	local DIR=$1
	echo_v 3 "Making directory: '$DIR'..."
	if [[ ! -e $DIR ]]; then
		mkdir $DIR
	elif [[ ! -d $DIR ]]; then
		echo "ERROR: '$DIR' already exists but is not a directory." 1>&2
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
	local EXT=$(GET_FILE_FROM_EXT $FULLFILE)
	EXT=${EXT^^}
	if [ $EXT == "HMM" ]; then 
		TYPE="HMM"
	elif [ $EXT == "FASTA" || $EXT == "FA" ]; then 
		TYPE="FASTA"
	elif [ $EXT == "HHM" ]; then 
		TYPE="HHM"
	elif [ $EXT == "MSA" ]; then 
		TYPE="MSA"
	elif [ $EXT == "MM_MSA" ]; then 
		TYPE="MM_MSA"
	elif [ $EXT == "MM_DB" ]; then 
		TYPE="MM_DB"
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

echo_v 3 "# mmore-functions.sh Loaded."
