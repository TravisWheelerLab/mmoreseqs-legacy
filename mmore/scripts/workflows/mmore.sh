#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore.sh	
#	AUTHOR:		David Rich
#	DESC: 		Run full MMORE pipeline.
###########################################################################

# PROGRAM:
# (1) verify proper arguments / filetypes
# (2) create mmseqs dbs
# (3) run mmseqs search
# (4) convert mmseqs output for mmore input
# (5) index mmore dbs
# (5) run mmore search

#set verbosity
VERBOSE="${VERBOSE:-3}"
if (( $VERBOSE > 3 )); then 
	VERB_ARG="-v"
fi

# ##### FUNCTIONS ###### #
# main script only contains the minimum functions needed to load helper function script
{
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
	echo_v 1 "Running Program: mmore.sh..."

	# load functions
	echo_v 3 "Importing 'mmore_functions.sh..."
	LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_functions.sh"
	# load tools 
	echo_v 3 "Importing 'mmore_get-tools.sh'..."
	LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_get-tools.sh"

	# variable preprocessing 
	{
		# process commandline variables
		{
			#verify proper number of variables
			NUM_ARGS=$#
			REQ_ARGS=4
			if (( $NUM_ARGS < $REQ_ARGS )); then 
				echo "ERROR: Illegal number of main args: ($NUM_ARGS of $REQ_ARGS)"
				echo "Usage: <target_mmore> <query_mmore> <target_mmseqs_p> <target_mmseqs_s> <query_mmseqs> | <target_type> <query_type> <target_mmseqs_type> <query_mmseqs_type>"
				echo "Accepted Filetypes: FASTA, HMM, HHM, MSA, MM_MSA, MMDB_S, MMDB_P"
				echo ""
				exit 1
			fi

			# main commandline args
			ARG_TARGET=$1
			ARG_QUERY=$2
			ARG_TARGET_MMSEQS_P=$3
			ARG_TARGET_MMSEQS_S=$4
			ARG_QUERY_MMSEQS=$5
			# optional commandline args
			ARG_TARGET_TYPE=$5
			ARG_QUERY_TYPE=$6
			ARG_TARGET_MMSEQS_P_TYPE=$7
			ARG_TARGET_MMSEQS_S_TYPE=$7
			ARG_QUERY_MMSEQS_TYPE=$8

			# print command line arguments
			echo_v 1 "./mmore.sh [0]$ARG_TARGET [1]$ARG_QUERY [2]$ARG_TARGET_MMSEQS [3]$ARG_QUERY_MMSEQS | [4]$ARG_TARGET_TYPE [5]$ARG_QUERY_TYPE [6]$ARG_TARGET_MMSEQS_TYPE [7]$ARG_QUERY_MMSEQS_TYPE"
		}

		# process environmental variables: if optional args not supplied by environment, falls back on these defaults
		{
			# main args:
			TARGET="${TARGET_MMORE:-$ARG_TARGET_MMORE}"
			QUERY="${QUERY_MMORE:-$ARG_QUERY_MMORE}"
			TARGET_MMSEQS_P="${TARGET_MMSEQS_P:-$ARG_TARGET_MMSEQS_P}"
			TARGET_MMSEQS_S="${TARGET_MMSEQS_S:-$ARG_TARGET_MMSEQS_S}"
			QUERY_MMSEQS="${QUERY_MMSEQS:-$ARG_QUERY_MMSEQS}"
			# main arg types:
			TARGET_TYPE="${TARGET_MMORE_TYPE:-$ARG_TARGET_MMORE_TYPE}"
			QUERY_TYPE="${QUERY_MMORE_TYPE:-$ARG_QUERY_MMORE_TYPE}"
			TARGET_MMSEQS_P_TYPE="${TARGET_MMSEQS_P_TYPE:-$ARG_TARGET_MMSEQS_P_TYPE}"
			TARGET_MMSEQS_S_TYPE="${TARGET_MMSEQS_S_TYPE:-$ARG_TARGET_MMSEQS_S_TYPE}"
			QUERY_MMSEQS_TYPE="${QUERY_MMSEQS_TYPE:-$ARG_QUERY_MMSEQS_TYPE}"
			# options :
			ROOT_DIR="${ROOT_DIR:-'/'}"
			PREP_DIR="${PREP_DIR:-""}"
			SCRIPT_DIR="${SCRIPT_DIR:-""}"
			TEMP_DIR="${TEMP_DIR:-""}"
			# pipeline options:
			RM_TEMP="${RM_TEMP:-0}" 
			DO_PREP="${DO_PREP:-0}"
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
			MAIN_FILES["TARGET_MMSEQS_P"]="$TARGET_MMSEQS_P" 
			MAIN_FILES["TARGET_MMSEQS_S"]="$TARGET_MMSEQS_S" 
			MAIN_FILES["QUERY_MMSEQS"]="$QUERY_MMSEQS"

			declare -A MAIN_TYPES
			MAIN_TYPES["TARGET_TYPE"]="$TARGET_TYPE"
			MAIN_TYPES["QUERY_TYPE"]="$QUERY_TYPE" 
			MAIN_TYPES["TARGET_MMSEQS_P_TYPE"]="$TARGET_MMSEQS_P_TYPE" 
			MAIN_TYPES["TARGET_MMSEQS_S_TYPE"]="$TARGET_MMSEQS_S_TYPE" 
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
			declare -A OUTPUT_FILES
		}

		# report variables
		{
			if (( $VERBOSE >= 3 )); then
				echo "#           TARGET: $TARGET"
				echo "#            QUERY: $QUERY"
				echo "#  TARGET_MMSEQS_P: $TARGET_MMSEQS_P"
				echo "#  TARGET_MMSEQS_S: $TARGET_MMSEQS_S"
				echo "#     QUERY_MMSEQS: $QUERY_MMSEQS"
				echo "#          RESULTS: $RESULTS"
				echo ""
				echo "#     FILTER_SCORE: $K_SCORE"
				echo "#      UNGAP_SCORE: $MIN_UNGAPPED_SCORE"
				echo "#        GAP_SCORE: $PVAL_CUTOFF => $EVAL_CUTOFF"
				echo ""
				echo "#            ALPHA: $ALPHA"
				echo "#             BETA: $BETA"
				echo "#            GAMMA: $GAMMA"
				echo "#     REPORT_SCORE: $PVALUE_REPORT => $EVALUE_REPORT"
			fi
		}
	}

	# === FILE PRE-PROCESSING === #
	{
		echo_v 0 "# (0/4) File formatting, verifying proper filetypes..."

		# build temporary folders
		{
			TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-mmore-XXXX)}"
			# top-level temporary directory
			TMP=${TEMP_DIR}/
			TMP_QUERY=${TMP}/query/
			TMP_TARGET=${TMP}/target/
			# temporary subdirectory for mmore
			TMP_MMORE=${TMP}/mmore/
			TMP_MMORE_DB=${TMP_MMORE}/db/
			TMP_MMORE_OUT=${TMP_MMORE}/out/
			# temporary subdirectory for mmseqs
			TMP_MMSEQS=${TMP}/mmseqs/
			TMP_MMSEQS_DB=${TMP_MMSEQS}/db/
			TMP_MMSEQS_WORKING=${TMP_MMSEQS}/working/
			TMP_MMSEQS_OUT=${TMP_MMSEQS}/out/

			# make tmp directories
			MAKE_DIR $TMP_ROOT
			MAKE_DIR $TMP_MMORE
			MAKE_DIR $TMP_MMORE_DB
			MAKE_DIR $TMP_MMORE_OUT
			MAKE_DIR $TMP_MMSEQS
			MAKE_DIR $TMP_MMSEQS_DB
			MAKE_DIR $TMP_MMSEQS_WORKING
			MAKE_DIR $TMP_MMSEQS_OUT

			# list of temp folders 
			TEMP_FOLDERS=""
			TEMP_FOLDERS+="${TMP}" 
			TEMP_FOLDERS+="${TMP_MMSEQS}" 
			TEMP_FOLDERS+="${TMP_MMORE}"
			TEMP_FOLDERS+="${TMP_MMSEQS_DB}"
		}

		# if type is not given, try to infer type based on file extension
		{
			if [ -z $TARGET_TYPE ]; then 
				TARGET_TYPE=$(INFER_FILETYPE )
				echo_v 3 "<TARGET_TYPE> inferred: $TARGET_TYPE"
			fi 
			if [ -z $QUERY_TYPE ]; then
				QUERY_TYPE=$(INFER_FILETYPE $QUERY)
				echo_v 3 "<QUERY_TYPE> inferred: $QUERY_TYPE"
			fi
			if [ -z $TARGET_MMSEQS_P_TYPE ]; then 
				TARGET_MMSEQS_P_TYPE=$(INFER_FILETYPE $TARGET_MMSEQS_P)
				echo_v 3 "<TARGET_MMSEQS_P_TYPE> inferred: $TARGET_MMSEQS_P_TYPE"
			fi 
			if [ -z $TARGET_MMSEQS_S_TYPE ]; then 
				TARGET_MMSEQS_TYPE_S=$(INFER_FILETYPE $TARGET_MMSEQS_S)
				echo_v 3 "<TARGET_MMSEQS_TYPE_S> inferred: $TARGET_MMSEQS_S_TYPE"
			fi 
			if [ -z $QUERY_MMSEQS_TYPE ]; then 
				QUERY_MMSEQS_TYPE=$(INFER_FILETYPE $QUERY_MMSEQS)
				echo_v 3 "<QUERY_MMSEQS_TYPE> inferred: $QUERY_MMSEQS_TYPE"
			fi 
		}

		# if prep is allowed, now converts target, query, and target_mmseqs to suitable filetypes
		if [ -z $DO_PREP ]; then
			echo "DO_PREP = TRUE"
			# TODO: should call back to 
			# LOAD_SOURCE ${SCRIPT_DIR}/mmore-prep.sh $TARGET $QUERY $TARGET_MMSEQS $TARGET_TYPE $QUERY_TYPE $TARGET_MMSEQS_TYPE
		fi

		# verify that input files are proper type
		{
			:
			# if [ $TARGET_TYPE != "HMM" ] 
			# then 
			# 	ERROR="TARGET_BAD_TYPE"
			# 	echo "ERROR: <TARGET> must be filetype: HMM (Given: $TARGET_TYPE). Can be converted from MSA or FASTA."
			# fi
			# if [ $QUERY_TYPE != "FASTA" ] 
			# then 
			# 	ERROR="QUERY_BAD_TYPE"
			# 	echo "ERROR: <QUERY> must be filetype: FASTA (Given: $QUERY_TYPE). Can be converted from HMM or MSA."
			# fi
			# if [ $TARGET_MMSEQS_P_TYPE != "HHM" ] || [ $TARGET_MMSEQS_P_TYPE != "FASTA" ] || [ $TARGET_MMSEQS_P_TYPE != "MMDB_P" ] || [ $TARGET_MMSEQS_P_TYPE != "MMDB_S" ]
			# then	
			# 	ERROR="TARGET_MMSEQS_BAD_TYPE"
			# 	echo "ERROR: <TARGET_MMSEQS> must be filetype: FASTA, or MMDB_S (Given: ${TARGET_MMSEQS_TYPE}). Can be converted from MSA or MM_MSA."
			# fi
			# if [ $QUERY_MMSEQS_TYPE != "FASTA" ] || [ $QUERY_MMSEQS_TYPE != "MMDB_S" ]
			# then
			# 	ERROR="QUERY_MMSEQS_BAD_TYPE"
			# 	echo "ERROR: <QUERY_MMSEQS> must be filetype: FASTA or MM_DB (Given: $QUERY_MMSEQS_TYPE). Can be converted from MSA or MM_MSA."
			# fi 
			# if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 1 )); 
			# then 
			# 	echo "HINT: To perform conversions, run again with 'mmore mmore -- prep 1' or run before with 'mmore prep'."
			# 	exit 1
			# fi

			TARGET_TYPE="HMM"
			QUERY_TYPE="FASTA"
			TARGET_MMSEQS_P_TYPE="MMDB_P"
			TARGET_MMSEQS_S_TYPE="MMDB_S"
			QUERY_MMSEQS_TYPE="MMDB_S"
		} 

		# Verify that input files exist
		{
			FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET)
			if [ -z $FILE_EXISTS ]; then 
				ERROR="TARGET_BAD_INPUT"
				echo_v 3 "ERROR: <TARGET> input '$TARGET' does not exist."
			fi 
			FILE_EXISTS=$(CHECK_FILE_EXISTS $QUERY)
			if [ -z $FILE_EXISTS ]; then
				ERROR="QUERY_BAD_INPUT"
				echo_v 3 "ERROR: <QUERY> input '$QUERY' does not exist."
			fi
			FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET_MMSEQS_S)
			if [ -z $FILE_EXISTS ]; then
				ERROR="TARGET_MMSEQS_BAD_INPUT"
				echo_v 3 "ERROR: <TARGET_MMSEQS_S> input '$TARGET_MMSEQS_S' does not exist."
			fi
			FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET_MMSEQS_P)
			if [ -z $FILE_EXISTS ]; then
				ERROR="TARGET_MMSEQS_BAD_INPUT"
				echo_v 3 "ERROR: <TARGET_MMSEQS_P> input '$TARGET_MMSEQS_P' does not exist."
			fi
			FILE_EXISTS=$(CHECK_FILE_EXISTS $QUERY_MMSEQS)
			if [ -z $FILE_EXISTS ]; then
				ERROR="QUERY_MMSEQS_BAD_INPUT"
				echo_v 3 "ERROR: <QUERY_MMSEQS> input '$QUERY_MMSEQS' does not exist."
			fi
			# if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 1 )); then 
			# 	echo "ERROR: Please check your filepaths."
			# 	exit 1
			# fi
		}

		# If overwriting is not permitted, check that destination files do not already exist
		{
			if (( $DO_OVERWRITE == 0 ))
			then
				OVERWRITE_FILES=""
				for FILE in $OUTPUT_FILES; do
					FILE_EXISTS=$(CHECK_FILE_EXISTS $FILE)
					if [ -z $FILE_EXISTS ]; then
						ERROR="OVERWRITE_FILE"
						OVERWRITE_FILES+="$FILE"
					fi 
				done

				if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 0 )); then 
					echo "ERROR: Output attempting to overwrite the following pre-existing file(s): ${OVERWRITE FILES[@]}"
					echo "Hint: To allow overwriting old files, use flag '--overwrite 1'"
					ERROR_EXIT 1
				fi
			fi 
		}

		# Determine type of search being run (Profile-to-Sequence or Sequence-to-Sequence)
		{
			# if mmseqs's main target database is  
			if [[ $TARGET_MMSEQS_P_TYPE == "MMDB_P" ]] || [[ $TARGET_MMSEQS_P_TYPE == "HHM" ]]
			then 
				SEARCH_TYPE="P2S"
			else
				SEARCH_TYPE="S2S"
			fi
		}

		# Soft link input files into temporary working folders 
		{
			:
		}
	}

	# ==== MMORE INDEX ==== #
	{
		echo "# (1/4) Create MMORE index of database..."
		# check if input files exist
		# if [ -f $TARGET_INDEX ]; then
		# 	# if so, move copy to temp folder
		# 	FILENAME=$(GET_FILE_FROM_PATH $TARGET_INDEX)
		# 	cp $TARGET_INDEX $TMP_MMORE/$FILENAME 
		# else
		# 	echo "build target index"
		# 	# if not, build index now 
		# 	$MMORE index $TARGET no-file --index $TARGET_INDEX no-file
		# fi

		# # check if input files exist
		# if [ -f $QUERY_INDEX ]; then
		# 	# if so, move copy to temp folder
		# 	FILENAME=$(GET_FILE_FROM_PATH $QUERY_INDEX)
		# 	cp $QUERY_INDEX $TMP_MMORE/$FILENAME 
		# else
		# 	echo "build query index"
		# 	# if not, build index now 
		# 	$MMORE index $TARGET no-file --index $TARGET_INDEX no-file
		# fi
		# TARGET_INDEX=$TMP_MMORE/target.idx 
		# QUERY_INDEX=$TMP_MMORE/query.idx

		# $MMORE index $TARGET $QUERY

		# mv $TARGET.idx $TARGET_INDEX 
		# mv $QUERY.idx $QUERY_INDEX

		# # capture database size
		# NUM_TARGETS=$( grep "#" -v $TARGET_INDEX | wc -l )
		# NUM_QUERIES=$( grep "#" -v $QUERY_INDEX | wc -l )
		# DB_SIZE=$((NUM_TARGETS))
		# echo "DB_SIZE: $DB_SIZE"
	}

	# ==== MMSEQS BUILD DATABASE === #
	{
		:
		# # create target profile database if input type is not already MMDB
		# {
		# 	if [ $TARGET_MMSEQS_P_TYPE != "MMDB_P" ] || [ $TARGET_MMSEQS_P_TYPE != "MMDB_S" ]
		# 	then
		# 		# database name
		# 		FILENAME=$(GET_FILE_NAME $TARGET_MMSEQS_P)
		# 		TARGET_MMSEQS_P_DB="${TMP_MMSEQS_DB}/${FILENAME}.TARGET_P.mmdb"
		# 		TARGET_MMSEQS_P_DB="${TMP_MMSEQS_DB}/target.p.mmdb"

		# 		# create database
		# 		$MMSEQS createdb 								\
		# 		$TARGET_MMSEQ_P $TARGET_MMSEQS_P_DB 	\

		# 		# check exitcode of process
		# 		EXIT_CODE=$?
		# 		if (( $EXIT_CODE != 0 && $DO_IGNORE_WARNINGS != 0 )); then
		# 			echo "ERROR: mmseqs createdb <TARGET_MMSEQS_PS> failed."
		# 			ERROR_EXIT $EXIT_CODE
		# 		fi

		# 		TARGET_MMSEQS_S=$TARGET_MMSEQS_P_DB
		# 		TARGET_MMSEQS_TYPE="MMDB_P"
		# 	fi
		# }

		# # create target sequence database if input type is not already MMDB
		# {
		# 	if [ $TARGET_MMSEQS_S_TYPE != "MMDB_S" ]
		# 	then
		# 		# database name
		# 		FILENAME=$(GET_FILE_NAME $TARGET_MMSEQS_S)
		# 		TARGET_MMSEQS_S_DB="${TMP_MMSEQS_DB}/${FILENAME}.TARGET_S.mmdb"
		# 		TARGET_MMSEQS_S_DB="${TMP_MMSEQS_DB}/target.s.mmdb"

		# 		# create database
		# 		$MMSEQS createdb 								\
		# 		$TARGET_MMSEQ_S $TARGET_MMSEQS_S_DB 	\

		# 		# check exitcode of process
		# 		EXIT_CODE=$?
		# 		if (( $EXIT_CODE != 0 && $DO_IGNORE_WARNINGS != 0 )); then
		# 			echo "ERROR: mmseqs createdb <TARGET_MMSEQS_S> failed."
		# 			ERROR_EXIT $EXIT_CODE
		# 		fi

		# 		TARGET_MMSEQS_S=$TARGET_MMSEQS_S_DB
		# 		TARGET_MMSEQS_S_TYPE="MMDB_S"
		# 	fi
		# }

		# # create query database if input type is not already MM_DB
		# {
		# 	if [ "QUERY_MMSEQS_TYPE" != "MMDB_S" ]
		# 	then
		# 		# database name
		# 		FILENAME=$(GET_FILE_NAME $QUERY_MMSEQS)
		# 		QUERY_MMSEQS_DB="${TMP_MMSEQS_DB}/${FILENAME}.QUERY.mmdb"
		# 		QUERY_MMSEQS_DB="${TMP_MMSEQS_DB}/query.mmdb"

		# 		# create database
		# 		$MMSEQS createdb 							\
		# 		$QUERY_MMSEQS $QUERY_MMSEQS_DB 	   \

		# 		# check exitcode of process
		# 		EXIT_CODE=$?
		# 		if (( $EXIT_CODE != 0 && $DO_IGNORE_WARNINGS != 0 )); then
		# 			echo "ERROR: mmseqs createdb <QUERY_MMSEQS> failed."
		# 			ERROR_EXIT $EXIT_CODE
		# 		fi

		# 		TARGET_MMSEQS=$TARGET_MMSEQS_DB
		# 		TARGET_MMSEQS_TYPE="MMDB_S"
		# 	fi
		# }

		MMSEQS_QUERY_DB=$QUERY_MMSEQS
		MMSEQS_TARGET_P_DB=$TARGET_MMSEQS_P 
		MMSEQS_TARGET_S_DB=$TARGET_MMSEQS_S
	}

	# ==== MMSEQS SEARCH ======= #
	{
		# Configure MMSEQS Options
		{
			PREFILTER_OPTS=""

			P2S_OPTIONS=""

			S2S_OPTIONS=""

			CONVERT_OPTIONS=""

			# --format-output STR       		CSV list of output columns from: query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid,taxid,taxname,taxlineage  [query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits]
			# -s FLOAT                     	Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [5.700]
			# -k INT                       	k-mer length (0: automatically set to optimum) [0]
			# --k-score INT                	k-mer threshold for generating similar k-mer lists [2147483647]
			# --max-rejected INT           	Maximum rejected alignments before alignment calculation for a query is stopped [2147483647]
			# --max-accept INT             	Maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
			# --local-tmp STR               	Path where some of the temporary files will be created []
			# --remove-tmp-files BOOL       	Delete temporary files [1]
			# --e-profile FLOAT             	Include sequences matches with < e-value thr. into the profile (>=0.0) [0.001]
		}

		# Run MMSEQS prefilter
		{
			MMSEQS_PREFILTER_DB=${TMP_MMSEQS_WORKING}/prefilter_aln.mmdb
			MMSEQS_PREFILTER_REPORT=${TMP}/mmseqs-prefilter.out 

			$MMSEQS prefilter 								\
			$MMSEQS_TARGET_P_DB $MMSEQS_QUERY_DB		\
			$MMSEQS_PREFILTER_DB 							\
			-v 			$VERBOSE 							\
			--k-score 	$K_SCORE 							\
			# > $MMSEQS_PREFILTER_REPORT 					\

			EXIT_CODE=$?
		}

		# Run MMSEQS Profile-to-Sequence alignment search
		{
			MMSEQS_P2S_DB=${TMP_MMSEQS_WORKING}/p2s_aln.mmdb 

			$MMSEQS align 										\
			$MMSEQS_TARGET_P_DB  $MMSEQS_QUERY_DB 		\
			$MMSEQS_PREFILTER_DB 							\
			$MMSEQS_P2S_DB 									\
			-v 			$VERBOSE 							\

			EXIT_CODE=$?
		}

		# Translate MMSEQS Profile IDs to MMSEQS Sequence IDs
		{
			echo "Running: translate_ids..."
			bash ${SCRIPT_DIR}/helpers/translate_ids.sh $MMSEQS_TARGET_P_DB $MMSEQS_TARGET_S_DB $MMSEQS_P2S_DB $TMP_MMSEQS_WORKING/dummy.mmdb 0
		}

		# Run MMSEQS Sequence-to-Sequence search
		{
			MMSEQS_S2S_DB=${TMP_MMSEQS_WORKING}/s2s_aln.mmdb

			$MMSEQS align 										\
			$MMSEQS_TARGET_S_DB $MMSEQS_QUERY_DB 		\
			$MMSEQS_PREFILTER_DB 							\
			$MMSEQS_S2S_DB 									\
			-v 			$VERBOSE 							\

			EXIT_CODE=$?
		}

		# Report results to m8 file
		{
			MMSEQS_M8="${TMP_MMSEQS_OUT}/mmseqs.results.m8"
			echo "mmseqs_results: $MMSEQS_M8"

			$MMSEQS convertalis 								\
			$MMSEQS_TARGET_S_DB $MMSEQS_QUERY_DB 		\
			$MMSEQS_S2S_DB 									\
			$MMSEQS_M8 	 										\

			EXIT_CODE=$?
		}
		
		# FORMAT_OUTPUT="qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
		# EVAL_CUTOFF=$(( PVAL_CUTOFF * DB_SIZE ))
		# EVAL_REPORT=$(( PVAL_REPORT * DB_SIZE ))

		# time $MMSEQS search 
		# 	$QUERY_MMSEQS  			$TARGET_MMSEQS 						\
		# 	$RESULT_MMSEQS  			$TMP_MMSEQS 							\
		# 	-k 							$KMER 									\
		# 	--min-ungapped-score 	$MIN_UNGAPPED_SCORE 					\
		# 	-e 							$EVAL_CUTOFF 							\
		# 	--k-score 					$K_SCORE 								\
		# 	--remove-tmp-files  		0											\
	}

	# # ==== CONVERT MMSEQS OUTPUT -> MMORE INPUT === #
	{
		:
	}

	# # ==== MMORE SEARCH ==== #
	{
		:
		echo "(7/7) Performing MMORE Search..."
		# cloud search options:
		# --alpha FLOAT 					alpha cloud tuning parameter
		# --beta INT 						beta cloud tuning parameter
		# --pipeline STR/INT 			Choose: [test,main,mmseqs,index]
		# --threshold FLOAT 				Score Reporting Threshold
		# --format-input STR 			CSV list of input format
		# --mmseqs-results STR 			Location of mmseqs results
		# --mmseqs-tid STR 				Location of mmseqs target id lookup table
		# --mmseqs-qid STR 				Location of mmseqs query id lookup table 

		TARGET_INDEX="${TMP_MMORE}/target.idx" 
		QUERY_INDEX="${TMP_MMORE}/query.idx"

		# Configure Options 
		{
			MMORE_OPTS=""

			# check if output are set
			if [ ! -z "$STDOUT" ]; then
				MMORE_OPTS="--stdout $STDOUT"
			fi
			if [ ! -z "$MMORE_M8OUT" ]; then
				M8_OUT=""
			else
				STDOUT="--m8-out $STDOUT"
			fi
		}

		# Run Search
		{
			echo "$MMORE mmore_main $TARGET $QUERY"
			MMORE_M8="${TMP_MMORE_OUT}/mmore.results.m8out"
			MMORE_MYTIME="${TMP_MMORE_OUT}/mmore.results.mytimeout"

			$MMORE mmore_main 											\
			$TARGET $QUERY 												\
			--mmseqs-m8 		$MMSEQS_M8	 							\
			--alpha 				$ALPHA 									\
			--beta 				$BETA 									\
			--gamma 				$GAMMA 									\
																				\
			--m8out 				$MMORE_M8 								\
			--mytimeout 		$MMORE_MYTIME 							\


		}

		mv $RESULT_CLOUD $RESULTS
	}

	# # === CLEAN UP === #
	# {
	# 	# move output files to their proper final destinations
	# 	{

	# 	}

	# 	# if temporary files flagged for removal, do it now
	# 	{
	# 		# TODO: remove temporary folders
	# 		# NOTE this should be done safely, removing files created by process, but not deleting folders unless they are empty
	# 		if (( $REMOVE_TEMP == 1 ))
	# 		then
	# 			for key in "${TEMP_FILES[@]}"; do
	# 				rm 
	# 		fi
	# 	}
	# }

}