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

}

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
		echo_v 1 "./mmore.sh [0]$ARG_TARGET [1]$ARG_QUERY [2]$ARG_TARGET_MMSEQS [3]$ARG_QUERY_MMSEQS | [4]$ARG_TARGET_TYPE [5]$ARG_QUERY_TYPE [6]$ARG_TARGET_MMSEQS_TYPE [7]$ARG_QUERY_MMSEQS_TYPE"
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
		TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-mmore-XXXX)}"
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

# === FILE PRE-PROCESSING === #
{
	echo_v 0 "# (0/4) File formatting, verifying proper filetypes..."

	# if type is not known, try to infer type based on file extension
	{
		if [ -z $TARGET_TYPE ]; then 
			TARGET_TYPE=$(INFER_TYPE )
			echo_v 3 "<TARGET_TYPE> inferred: $TARGET_TYPE"
		fi 
		if [ -z $QUERY_TYPE ]; then
			QUERY_TYPE=$(INFER_TYPE $QUERY)
			echo_v 3 "<QUERY_TYPE> inferred: $QUERY_TYPE"
		fi
		if [ -z $TARGET_MMSEQS_TYPE ]; then 
			TARGET_MMSEQS_TYPE=$(INFER_TYPE $TARGET_MMSEQS)
			echo_v 3 "<TARGET_MMSEQS_TYPE> inferred: $TARGET_MMSEQS_TYPE"
		fi 
		if [ -z $QUERY_MMSEQS_TYPE ]; then 
			QUERY_MMSEQS_TYPE=$(INFER_TYPE $QUERY_MMSEQS)
			echo_v 3 "<QUERY_MMSEQS_TYPE> inferred: $QUERY_MMSEQS_TYPE"
		fi 
	}

	# if prep is allowed, now converts target, query, and target_mmseqs to suitable filetypes
	if [ -z $DO_PREP ]; then
		# TODO: should call back to 
		LOAD_SOURCE ${SCRIPT_DIR}/mmore-prep.sh $TARGET $QUERY $TARGET_MMSEQS $TARGET_TYPE $QUERY_TYPE $TARGET_MMSEQS_TYPE
	fi

	# verify that target, query, and target_mmseqs are suitable filetypes
	{
		if [ $TARGET_TYPE != "HMM" ] 
		then 
			ERROR="TARGET_BAD_TYPE"
			echo "ERROR: <TARGET> must be filetype HMM (Given: $TARGET_TYPE). Can be converted from MSA or FASTA."
		fi
		if [ $QUERY_TYPE != "FASTA" ] 
		then 
			ERROR="QUERY_BAD_TYPE"
			echo "ERROR: <QUERY> must be filetype FASTA (Given: $QUERY_TYPE). Can be converted from HMM or MSA."
		fi
		if [ $TARGET_MMSEQS_TYPE != "HHM" ] || [ $TARGET_MMSEQS_TYPE != "FASTA" ] || [ $TARGET_MMSEQS_TYPE != "MM_DB" ]
		then	
			ERROR="TARGET_MMSEQS_BAD_TYPE"
			echo "ERROR: <TARGET_MMSEQS> must be filetype HHM, FASTA, or MM_DB (Given: ${TARGET_MMSEQS_TYPE}). Can be converted from MSA or MM_MSA."
		fi 
		if [ $QUERY_MMSEQS_TYPE != "FASTA" ] || [ $QUERY_MMSEQS_TYPE != "MM_DB" ]
		then
			ERROR="QUERY_MMSEQS_BAD_TYPE"
			echo "ERROR: <QUERY_MMSEQS> must be filetype FASTA or MM_DB (Given: $QUERY_MMSEQS_TYPE). Can be converted from MSA or MM_MSA."
		fi 
		if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 1 )); 
		then 
			echo "HINT: To perform conversions, run again with 'mmore mmore -- prep 1' or run before with 'mmore prep'."
			exit 1
		fi
	}

	# verify that main files exist
	{
		FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET)
		if [ ! -z $FILE_EXISTS ]; then 
			ERROR="TARGET_BAD_INPUT"
			echo_v 3 "ERROR: <TARGET> input '$TARGET' does not exist."
		fi 
		FILE_EXISTS=$(CHECK_FILE_EXISTS $QUERY)
		if [ -z $FILE_EXISTS ]; then
			ERROR="QUERY_BAD_INPUT"
			echo_v 3 "ERROR: <QUERY> input '$QUERY' does not exist."
		fi
		FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET_MMSEQS)
		if [ -z $FILE_EXISTS ]; then
			ERROR="TARGET_MMSEQS_BAD_INPUT"
			echo_v 3 "ERROR: <TARGET_MMSEQS> input '$TARGET_MMSEQS' does not exist."
		fi
		FILE_EXISTS=$(CHECK_FILE_EXISTS $QUERY_MMSEQS)
		if [ -z $FILE_EXISTS ]; then
			ERROR="QUERY_MMSEQS_BAD_INPUT"
			echo_v 3 "ERROR: <QUERY_MMSEQS> input '$QUERY_MMSEQS' does not exist."
		fi
		if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 1 )); then 
			echo "ERROR: Please check your filepaths."
			exit 1
		fi
	}

	# TODO: if overwriting is not permitted, check that destination files do not already exist
	if (( $DO_OVERWRITE == 0 ))
	{
		if [ -z $ERROR ] && (( $DO_IGNORE_WARNINGS != 0 )); then 
			echo "Hint: to allow overwriting old files, use flag '--overwrite 1'"
			ERROR_EXIT 1
		fi
	}
}

# ==== MMORE INDEX ==== #
{
	echo "# (1/4) Create MMORE index of database..."
	# check if input files exist
	if [ -f $TARGET_INDEX ]; then
		# if so, move copy to temp folder
		FILENAME=$(GET_FILE_FROM_PATH $TARGET_INDEX)
		cp $TARGET_INDEX $TMP_MMORE/$FILENAME 
	else
		# if not, build index now 
	fi

	# check if input files exist
	if [ -f $QUERY_INDEX ]; then
		# if so, move copy to temp folder
		FILENAME=$(GET_FILE_FROM_PATH $QUERY_INDEX)
		cp $QUERY_INDEX $TMP_MMORE/$FILENAME 
	else
		# if not, build index now 
	fi

	# time $CLOUD index $TARGET $QUERY
	# time mv $TARGET.idx $TMP_MMORE/target.idx 
	# time mv $QUERY.idx $TMP_MMORE/query.idx

	# capture database size
	NUM_TARGETS=$( grep "#" -v $TMP_CLOUD/target.idx | wc -l )
	NUM_QUERIES=$( grep "#" -v $TMP_CLOUD/query.idx | wc -l )
	DB_SIZE=$((NUM_TARGETS))
}

# ==== MMSEQS BUILD DATABASE === #
{
	# create target database if input type is not already MM_DB
	{
		if [ $TARGET_MMSEQS_TYPE != "MM_DB" ]
		then
			FILENAME=$(GET_FILE_NAME $TARGET_MMSEQS)
			TARGET_MMSEQS_DB="${TMP_MMSEQS_DB}/${FILENAME}.TARGET.mm_db"

			# create database
			$MMSEQS createdb 							\
			$TARGET_MMSEQS $TARGET_MMSEQS_DB 	\

			# check exitcode of process
			EXIT_CODE=$?
			if (( $EXIT_CODE != 0 )) && (( $DO_IGNORE_WARNINGS != 0 )) then;
				echo "ERROR: mmseqs createdb <TARGET_MMSEQS> failed."
				ERROR_EXIT $EXIT_CODE
			fi

			TARGET_MMSEQS=$TARGET_MMSEQS_MMDB
			TARGET_MMSEQS_TYPE="MM_DB"
		fi
	}

	# create query database if input type is not already MM_DB
	{
		if [ "QUERY_MMSEQS_TYPE" != "MM_DB" ]
		then
			FILENAME=$(GET_FILE_NAME $QUERY_MMSEQS)
			QUERY_MMSEQS_DB="${TMP_MMSEQS_DB}/${FILENAME}.QUERY.MM_DB"

			# create database
			$MMSEQS createdb 							\
			$QUERY_MMSEQS $QUERY_MMSEQS_DB 	   

			# check exitcode of process
			EXIT_CODE=$?
			if (( EXIT_CODE != 0 )) then;
				echo "ERROR: mmseqs createdb <QUERY_MMSEQS> failed."
				ERROR_EXIT $EXIT_CODE
			fi

			TARGET_MMSEQS=$TARGET_MMSEQS_MMDB
			TARGET_MMSEQS_TYPE="MM_DB"
		fi
	}
	
}

# # ==== MMSEQS SEARCH ======= #
# {
# 	# Run MMSEQS prefilter
# 	{
# 		mmseqs prefilter 								\
# 		$MMSEQS_TARGET_DB $MMSEQS_QUERY_DB		\
# 		$MMSEQS_PREFILTER_DB 						\
# 		-k 	$K_SCORE 								\
# 		-v 	$VERBOSE 								\

# 		EXIT_CODE=$?
# 	}

# 	# Run MMSEQS Profile-to-Sequence search
# 	{

# 	}

# 	# Run MMSEQS Sequence-to-Sequence search
# 	{

# 	}

# 	# Map MMSEQS Profile IDs to MMSEQS Sequence IDs
# 	{

# 	}
	
# 	# --format-output STR       		CSV list of output columns from: query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid,taxid,taxname,taxlineage  [query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits]
# 	# -s FLOAT                     	Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [5.700]
# 	# -k INT                       	k-mer length (0: automatically set to optimum) [0]
# 	# --k-score INT                	k-mer threshold for generating similar k-mer lists [2147483647]
# 	# --max-rejected INT           	Maximum rejected alignments before alignment calculation for a query is stopped [2147483647]
# 	# --max-accept INT             	Maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
# 	# --local-tmp STR               	Path where some of the temporary files will be created []
# 	# --remove-tmp-files BOOL       	Delete temporary files [1]
# 	# --e-profile FLOAT             	Include sequences matches with < e-value thr. into the profile (>=0.0) [0.001]

# 	FORMAT_OUTPUT="qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
# 	EVAL_CUTOFF=$(( PVAL_CUTOFF * DB_SIZE ))
# 	EVAL_REPORT=$(( PVAL_REPORT * DB_SIZE ))

# 	time $MMSEQS search 
# 		$QUERY_MMSEQS  			$TARGET_MMSEQS 						\
# 		$RESULT_MMSEQS  			$TMP_MMSEQS 							\
# 		-k 							$KMER 									\
# 		--min-ungapped-score 	$MIN_UNGAPPED_SCORE 					\
# 		-e 							$EVAL_CUTOFF 							\
# 		--k-score 					$K_SCORE 								\
# 		--remove-tmp-files  		0											\
# }

# # ==== CONVERT MMSEQS OUTPUT -> MMORE INPUT === #

# # ==== MMORE SEARCH ==== #
# {
# 	echo "(7/7) Performing MMORE Search..."
# 	# cloud search options:
# 	# --alpha FLOAT 					alpha cloud tuning parameter
# 	# --beta INT 						beta cloud tuning parameter
# 	# --pipeline STR/INT 			Choose: [test,main,mmseqs,index]
# 	# --threshold FLOAT 				Score Reporting Threshold
# 	# --format-input STR 			CSV list of input format
# 	# --mmseqs-results STR 			Location of mmseqs results
# 	# --mmseqs-tid STR 				Location of mmseqs target id lookup table
# 	# --mmseqs-qid STR 				Location of mmseqs query id lookup table 

# 	TARGET_INDEX="$TMP_CLOUD/target.idx" 
# 	QUERY_INDEX="$TMP_CLOUD/query.idx"

# 	# create flags: for each option, if value is set, then generate flag
# 	{

# 	}

# 	# check if output are set
# 	if [ ! -z "$STDOUT" ]; then
# 		STDOUT=""
# 	else
# 		STDOUT="--stdout $STDOUT"
# 	fi
# 	if [ ! -z "$M8_OUT" ]; then
# 		M8_OUT=""
# 	else
# 		$STDOUT="--m8-out $STDOUT"
# 	fi
# 	if [ ! -z "$STDOUT" ]; then
# 		$STDOUT=""
# 	else
# 		$STDOUT="--stdout $STDOUT"
# 	fi

# 	time $MMORE mmore_main 											\
# 		$TARGET_PATH $QUERY_PATH 									\
# 		--mmseqs-m8 		$MMSEQS_M8OUT							\
# 		--alpha 				$ALPHA 									\
# 		--beta 				$BETA 									\
# 		--gamma 				$GAMMA 									\

# 	mv $RESULT_CLOUD $RESULTS
# }

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
