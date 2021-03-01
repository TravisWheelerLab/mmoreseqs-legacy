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
# (3a) if <target_mmseqs_p> is .msa file, convert to .mm_msa file.
# (3b) if <target_mmseqs_p> is .mm_msa file, convert to .hhm file.
# (3c) if <target_mmseqs_p> is .hhm file, convert to .mm_db file.
# (3d) if <target_mmseqs_p> is .fasta file, convert to .mm_db file.
# (3e) verify that final <target_mmseqs> file is .mm_db file.
# 
# (4a) if <query_mmseqs> is .msa file, convert to .fasta (consensus) file.
# (4b) if <query_mmseqs> is .mm_msa file, convert to .fasta (consensus) file (?).
# (4c) if <query_mmseqs> is .hhm file, convert to .fasta (consensus) file (?).
# (4d) if <query_mmseqs> is .fasta file, convert to .mm_db file.
# (4e) verify that final <query_mmseqs> file is .mm_db file.
#

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
{
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
			REQ_ARGS=3
			if (( $NUM_ARGS < $REQ_ARGS )); then 
				echo "ERROR: Illegal number of main args: ($NUM_ARGS of $REQ_ARGS)"
				echo "Usage: <target> <query> | <prep_folder> <target_type> <query_type>"
				echo "Accepted Filetypes: FASTA, MSA"
				exit 1
			fi

			# main commandline args
			ARG_TARGET=$1
			ARG_QUERY=$2
			# optional args
			ARG_TEMP_DIR=$3
			ARG_TARGET_TYPE_S=$4
			ARG_QUERY_TYPE=$5

			# print command line arguments
			echo_v 1 "./mmore.sh [0]$ARG_TARGET [1]$ARG_QUERY [2]$ARG_TARGET_MMSEQS [3]$ARG_QUERY_MMSEQS | [4]$ARG_TARGET_TYPE [5]$ARG_QUERY_TYPE [6]$ARG_TARGET_MMSEQS_TYPE [7]$ARG_QUERY_MMSEQS_TYPE"
		}

		# process environmental variables: if optional args not supplied by environment, falls back on these defaults
		{
			# main args:
			TARGET="${TARGET:-$ARG_TARGET}"
			QUERY="${QUERY:-$ARG_QUERY}"
			TEMP_DIR="${TEMP_DIR:-$ARG_TEMP_DIR}"
			TARGET_TYPE=
			# mmore args:
			TARGET_MMORE="${TARGET_MMORE:-$ARG_TARGET_MMORE}"
			QUERY_MMORE="${QUERY_MMORE:-$ARG_QUERY_MMORE}"
			TARGET_MMORE_TYPE="${TARGET_MMORE_TYPE:-$ARG_TARGET_MMORE_TYPE}"
			QUERY_MMORE_TYPE="${QUERY_MMORE_TYPE:-$ARG_QUERY_MMORE_TYPE}"
			# mmseqs args:
			TARGET_MMSEQS_P="${TARGET_MMSEQS_P:-$ARG_TARGET_MMSEQS_P}"
			TARGET_MMSEQS_S="${TARGET_MMSEQS_S:-$ARG_TARGET_MMSEQS_S}"
			QUERY_MMSEQS="${QUERY_MMSEQS:-$ARG_QUERY_MMSEQS}"
			TARGET_MMSEQS_P_TYPE="${TARGET_MMSEQS_P_TYPE:-$ARG_TARGET_MMSEQS_P_TYPE}"
			TARGET_MMSEQS_S_TYPE="${TARGET_MMSEQS_S_TYPE:-$ARG_TARGET_MMSEQS_S_TYPE}"
			QUERY_MMSEQS_TYPE="${QUERY_MMSEQS_TYPE:-$ARG_QUERY_MMSEQS_TYPE}"
			# directories:
			ROOT_DIR="${ROOT_DIR:-'/'}"
			SCRIPT_DIR="${SCRIPT_DIR:-""}"
			TEMP_DIR="${TEMP_DIR:-""}"
			# pipeline options:
			RM_TEMP="${RM_TEMP:-0}" 
			DO_PREP="${DO_PREP:-0}"
			DO_COPY="${DO_COPY:-1}"
			DO_OVERWRITE="${DO_OVERWRITE:-1}"
			DO_IGNORE_WARNINGS="${DO_IGNORE_WARNINGS:-1}"
			NUM_THREADS="${NUM_THREADS:-1}"
			# mmseqs options:
			KMER="${KMER:-7}"
			K_SCORE="${K_SCORE:-80}"
			MIN_UNGAPPED_SCORE="${UNGAPPEDVIT_THRESH:-15}"
			PVAL_CUTOFF="${GAPPEDVIT_THRESH:-0.001}"
			EVAL_CUTOFF="${EVAL_CUTOFF:-200}"
			SENSITIVITY="${SENSITIVITY:-7.5}"
			# mmore options:
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
			TMP="${TEMP_DIR}/"
			TMP_DB="${TMP}/db/"
			# mmore folders
			TMP_MMORE="${TMP}/mmore/"
			TMP_MMORE_DB="${TMP_MMORE}/db/"
			# mmseqs folders
			TMP_MMSEQS="${TMP}/mmseqs/"
			TMP_MMSEQS_DB="${TMP_MMSEQS}/db/"

			# make tmp directories
			MAKE_DIR $TMP
			MAKE_DIR $TMP_DB 
			MAKE_DIR $TMP_MMORE
			MAKE_DIR $TMP_MMORE_DB
			MAKE_DIR $TMP_MMSEQS
			MAKE_DIR $TMP_MMSEQS_DB

			# list of temp folders 
			TEMP_FOLDERS=""
			TEMP_FOLDERS+="${TMP}" 
			TEMP_FOLDERS+="${TMP_DB}"
			TEMP_FOLDERS+="${TMP_MMORE}"
			TEMP_FOLDERS+="${TMP_MMORE_DB}"
			TEMP_FOLDERS+="${TMP_MMSEQS}" 
			TEMP_FOLDERS+="${TMP_MMSEQS_DB}"
		}
	}

	# === FILE FORMATTING === #
	{
		echo "# (0/4) File formatting and building..."

		# COPY main files
		{
			# import using copy or symbolic links
			{
				LINK="ln -s"
				COPY="cp"

				if [ "$DO_COPY" == "0" ]
				then
					IMPORT="$LINK"
				else
					IMPORT="$COPY"
				fi
			}

			# verify valid target formats 
			if [ "$TARGET_TYPE" != "MSA" ] && [ "$TARGET_TYPE" != "FASTA" ]
			then 
				echo "ERROR: TARGET_TYPE is not a valid file type."
				exit 1
			fi 
			
			# import target file 
			TARGET_IN_TYPE="$TARGET_TYPE"
			if [ "$TARGET_TYPE" == "MSA" ]; then
				TARGET_MSA="${TMP_DB}/target.msa"
				$IMPORT "${TARGET}" "${TARGET_MSA}"
			elif [ "$TARGET_TYPE" == "FASTA" ]
				TARGET_FASTA="${TMP_DB}/target.fasta"
				$IMPORT "${TARGET}" "${TARGET_FASTA}"
			fi 

			# verify valid query formats 
			if [ "$QUERY_TYPE" != "MSA" ] && [ "$QUERY_TYPE" != "FASTA" ]
			then 
				echo "ERROR: QUERY_TYPE is not a valid file type."
				exit 1
			fi 

			# import query file
			QUERY_IN_TYPE="$QUERY_TYPE"
			if [ "$QUERY_TYPE" == "MSA" ]; then
				QUERY_MSA="${TMP_DB}/query.msa"
				$IMPORT "${QUERY}" "${QUERY_MSA}"
			elif [ "$QUERY_TYPE" == "FASTA" ]
				QUERY_FASTA="${TMP_DB}/query.fasta"
				$IMPORT "${QUERY}" "${QUERY_FASTA}"
			fi
		}

		# QUERY file formatting
		{
			# If query is a msa file, then convert it to hmm
			if [ "$QUERY_IN_TYPE" == "MSA" ]
			then
				QUERY_HMM="${TMP_DB}/query.hmm"
				$HMMBUILD "$QUERY_HMM" "$QUERY_MSA"

				QUERY_FASTA="${TMP_DB}/query.fasta"
				$HMMEMIT -c "$QUERY_HMM" > "$QUERY_FASTA"	
			fi

			# if target mmseqs file is a hhm, then we are ready.
			if [ "$QUERY_IN_TYPE" == "FASTA" ]
			then
				QUERY_HMM="${TMP_DB}/query.fasta"
				${SCRIPT_DIR}/helpers/convert_fasta-to-hmm.sh  
			fi

			# Final check if file is a valid file format
			if [ "$QUERY_TYPE" != "FASTA" ]
			then
				echo "ERROR: QUERY is an unsupported filetype of '$QUERY_TYPE'."
				echo "Valid QUERY Types: FASTA, MSA."
				exit 1
			fi
		}

		# QUERY_MMORE file formatting
		{
			
		}

		# QUERY_MMSEQS file formatting 
		{

		}

		# TARGET file formatting
		{
			TARGET_TYPE="$TARGET_IN_TYPE"

			# If target is a msa file 
			if [ "$TARGET_TYPE" == "MSA" ]
			then
				TARGET_HMM="${TMP_DB}/target.hmm"
				$HMMBUILD "$TARGET_HMM" "$TARGET_MSA"

				TARGET_TYPE="HMM"
			fi

			# If target is a hmm file, then convert to consensus fasta file
			if [ "$TARGET_TYPE" == "HMM" ]
			then
				TARGET_FASTA="${TMP_DB}/target.fasta"
				$HMMEMIT -c "$TARGET_HMM" > "$TARGET_FASTA"	

				TARGET_TYPE="FASTA"
			fi

			# if target mmseqs file is a hhm, then we are ready.
			if [ "$TARGET_TYPE" == "FASTA" ]
			then
				:
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

		# TARGET_MMORE file formatting 
		{

		}
	}
}
