#!/usr/bin/bash
###########################################################################
#	 - FILE:  mmore-search.sh	
#	 - DESC:  Run full MMORE pipeline.
# 	        This is the old version.
###########################################################################

# PROGRAM:
# (1) verify proper arguments / filetypes
# (2) create mmseqs dbs
# (3) run mmseqs search
# 	(3A) run mmseqs prefilter (P2S)
# 	(3B) run mmseqs align (P2S)
# 	(3C) run mmseqs align (S2S)
# 	(3D) sort and splice P2S and S2S results
# (4) convert mmseqs output for mmore input
# (5) index mmore dbs
# (5) run mmore search

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
	echo_v 2 "Running Program: mmore.sh..."
	PROGRAM+="> mmore-search.sh >"

	# load script dependencies 
	{
		# load functions
		echo_v 3 "Importing 'mmore_functions.sh..."
		LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_functions.sh"
		# load tools 
		# echo_v 3 "Importing 'mmore_get-tools.sh'..."
		# LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_get-tools.sh"
	}

	local BEG_TOTAL=$(GET_TIME)

	# variable preprocessing 
	{
		# process commandline variables
		{
			#verify proper number of variables
			NUM_ARGS=$#
			REQ_ARGS=5
			if (( $NUM_ARGS < $REQ_ARGS )); then 
				echo "ERROR: Illegal number of main args: ($NUM_ARGS of $REQ_ARGS)"
				echo "Usage: <i:target_mmore_p> <i:target_mmore_s> <i:query_mmore> <i:target_mmseqs_p> <i:target_mmseqs_s> <i:query_mmseqs> | <i:prep_dir>"
				echo "./mmore.sh [0]$ARG_TARGET [1]$ARG_QUERY [2]$ARG_TARGET_MMSEQS_P [3]$ARG_TARGET_MMSEQS_P [4]$ARG_QUERY_MMSEQS | [5]$PREP_FILE"
				echo "Accepted Filetypes: FASTA, HMM, HHM, MSA, MM_MSA, MMDB_S, MMDB_P"
				echo ""
				exit 1
			fi

			# main commandline args
			ARG_TARGET_MMORE=$1
			ARG_QUERY_MMORE=$2
			ARG_TARGET_MMSEQS_P=$3
			ARG_TARGET_MMSEQS_S=$4
			ARG_QUERY_MMSEQS=$5
			# optional commandline args
			ARG_PREP_DIR=$7
		}

		# process environmental variables: if optional args not supplied by environment, falls back on these defaults
		{
			# main args:
			TARGET_MMORE_P="${TARGET_MMORE_P:-$ARG_TARGET_MMORE_P}"
			TARGET_MMORE_S="${TARGET_MMORE_S:-$ARG_TARGET_MMORE_S}"
			QUERY_MMORE="${QUERY_MMORE:-$ARG_QUERY_MMORE}"
			TARGET_MMSEQS_P="${TARGET_MMSEQS_P:-$ARG_TARGET_MMSEQS_P}"
			TARGET_MMSEQS_S="${TARGET_MMSEQS_S:-$ARG_TARGET_MMSEQS_S}"
			QUERY_MMSEQS="${QUERY_MMSEQS:-$ARG_QUERY_MMSEQS}"
			# main arg types:
			PREP_DIR="${PREP_DIR:-$ARG_PREP_DIR}"
			TARGET_MMORE_P_TYPE="${TARGET_MMORE_P_TYPE:-$ARG_TARGET_MMORE_P_TYPE}"
			TARGET_MMORE_S_TYPE="${TARGET_MMORE_S_TYPE:-$ARG_TARGET_MMORE_S_TYPE}"
			QUERY_MMORE_TYPE="${QUERY_MMORE_TYPE:-$ARG_QUERY_MMORE_TYPE}"
			TARGET_MMSEQS_P_TYPE="${TARGET_MMSEQS_P_TYPE:-$ARG_TARGET_MMSEQS_P_TYPE}"
			TARGET_MMSEQS_S_TYPE="${TARGET_MMSEQS_S_TYPE:-$ARG_TARGET_MMSEQS_S_TYPE}"
			QUERY_MMSEQS_TYPE="${QUERY_MMSEQS_TYPE:-$ARG_QUERY_MMSEQS_TYPE}"
			
			# default args:
			SET_ENV_ARG_DEFAULTS
			SET_TOOLS
		}

		# list of all main arguments
		{
			# main files array
			declare -A MAIN_FILES
			MAIN_FILES["TARGET_MMORE"]="$TARGET_MMORE"
			MAIN_FILES["QUERY_MMORE"]="$QUERY_MMORE" 
			MAIN_FILES["TARGET_MMSEQS_P"]="$TARGET_MMSEQS_P" 
			MAIN_FILES["TARGET_MMSEQS_S"]="$TARGET_MMSEQS_S" 
			MAIN_FILES["QUERY_MMSEQS"]="$QUERY_MMSEQS"
			# main filetypes array
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
			if (( $VERBOSE >= 3 )); 
			then
				echo "# ========== MMORESEQS: SEARCH ================="
				echo "#     TARGET_MMORE: $TARGET_MMORE"
				echo "#      QUERY_MMORE: $QUERY_MMORE"
				echo "#  TARGET_MMSEQS_P: $TARGET_MMSEQS_P"
				echo "#  TARGET_MMSEQS_S: $TARGET_MMSEQS_S"
				echo "#     QUERY_MMSEQS: $QUERY_MMSEQS"
				echo "#"
				echo "#         TEMP_DIR: $TEMP_DIR"
				echo "#         PREP_DIR: $PREP_DIR"
				echo "#"
				echo "#    MMSEQS_KSCORE: $MMSEQS_KSCORE"
				echo "#  MMSEQS_UNGAPPED: $MMSEQS_UNGAPPED"
				echo "#      MMSEQS_PVAL: $MMSEQS_PVAL"
				echo "#"
				echo "#      MMORE_ALPHA: $MMORE_ALPHA"
				echo "#       MMORE_BETA: $MMORE_BETA"
				echo "#      MMORE_GAMMA: $MMORE_GAMMA"
				echo "#     MMORE_FILTER: $MMORE_DO_FILTER"
				echo "#    MMORE_VITERBI: $MMORE_VITERBI"
				echo "#      MMORE_CLOUD: $MMORE_CLOUD"
				echo "#   MMORE_BOUNDFWD: $MMORE_BOUNDFWD"
				echo "# =========================================="
			fi
		}
	}

	# === FILE PRE-PROCESSING === #
	{
		# build temporary folders
		{
			# assign temporary file if none given
			TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-mmore-XXXX)}"
			# top-level temporary directory
			TMP="${TEMP_DIR}/"
			TMP_QUERY="${TMP}/query/"
			TMP_TARGET="${TMP}/target/"
			# temporary subdirectory for mmore
			TMP_MMORE="${TMP}/mmore/"
			TMP_MMORE_DB="${TMP_MMORE}/db/"
			TMP_MMORE_OUT="${TMP_MMORE}/out/"
			# temporary subdirectory for mmseqs
			TMP_MMSEQS="${TMP}/mmseqs/"
			TMP_MMSEQS_DB="${TMP_MMSEQS}/db/"
			TMP_MMSEQS_WORKING="${TMP_MMSEQS}/working/"
			TMP_MMSEQS_OUT="${TMP_MMSEQS}/out/"

			# make tmp directories
			MAKE_DIR "$TMP"
			MAKE_DIR "$TMP_MMORE"
			MAKE_DIR "$TMP_MMORE_DB"
			MAKE_DIR "$TMP_MMORE_OUT"
			MAKE_DIR "$TMP_MMSEQS"
			MAKE_DIR "$TMP_MMSEQS_DB"
			MAKE_DIR "$TMP_MMSEQS_WORKING"
			MAKE_DIR "$TMP_MMSEQS_OUT"

			# list of temp folders 
			TEMP_FOLDERS=""
			TEMP_FOLDERS+="${TMP}" 
			TEMP_FOLDERS+="${TMP_MMSEQS}" 
			TEMP_FOLDERS+="${TMP_MMORE}"
			TEMP_FOLDERS+="${TMP_MMSEQS_DB}"
		}

		# if prep is allowed, now converts target, query, and target_mmseqs to suitable filetypes
		if [ -z $DO_PREP ]; then
			echo "DO_PREP = TRUE"
			# TODO: should call back to 
			# LOAD_SOURCE ${SCRIPT_DIR}/mmore-prep.sh $TARGET $QUERY $TARGET_MMSEQS $TARGET_TYPE $QUERY_TYPE $TARGET_MMSEQS_TYPE
		fi

		# verify that input files are proper type (TODO)
		{
			TARGET_TYPE="HMM"
			QUERY_TYPE="FASTA"
			TARGET_MMSEQS_P_TYPE="MMDB_P"
			TARGET_MMSEQS_S_TYPE="MMDB_S"
			QUERY_MMSEQS_TYPE="MMDB_S"
		} 

		# Verify that input files exist
		{
			FILE_EXISTS=$(CHECK_FILE_EXISTS $TARGET_MMORE)
			if [ -z $FILE_EXISTS ]; then 
				ERROR="TARGET_BAD_INPUT"
				echo_v 3 "ERROR: <TARGET> input '$TARGET_MMORE' does not exist."
			fi 
			FILE_EXISTS=$(CHECK_FILE_EXISTS $QUERY_MMORE)
			if [ -z $FILE_EXISTS ]; then
				ERROR="QUERY_BAD_INPUT"
				echo_v 3 "ERROR: <QUERY> input '$QUERY_MMORE' does not exist."
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
			# if [ -z $ERROR ] && (( $DO_IGNORE_ERRORS != 1 )); then 
			# 	echo "ERROR: Please check your filepaths."
			# 	exit 1
			# fi
		}

		# Interrim output files 
		{
			MMSEQS_PREFILTER="${TMP_MMSEQS_WORKING}/prefilter_aln.mmdb"
			MMSEQS_PREFILTER_STDOUT="${TMP}/mmseqs-prefilter.stdout"

			MMSEQS_P2S="${TMP_MMSEQS_WORKING}/p2s_aln.mmdb" 
			MMSEQS_P2S_STDOUT="${TMP_MMSEQS_OUT}/mmseqs-p2s.stdout" 

			MMSEQS_DUMMY="${TMP_MMSEQS_WORKING}/dummy.mmdb"
			MMSEQS_DUMMY_STDOUT="${TMP_MMSEQS_OUT}/dummy.stdout"

			MMSEQS_S2S="${TMP_MMSEQS_WORKING}/s2s_aln.mmdb"
			MMSEQS_S2S_STDOUT="${TMP_MMSEQS_OUT}/mmseqs-s2s.stdout"

			MMSEQS_P2S_M8="${TMP_MMSEQS_OUT}/mmseqs.p2s.mm_m8"
			MMSEQS_S2S_M8="${TMP_MMSEQS_OUT}/mmseqs.s2s.mm_m8"
			MMSEQS_M8="${TMP_MMSEQS_OUT}/mmseqs.mm_m8"
			MMSEQS_CONVERT_STDOUT="${TMP_MMSEQS_OUT}/mmseqs-convertalis.stdout"

			MMORE_STDOUT="${TMP_MMORE_OUT}/mmore.results.stdout"
			MMORE_STDERR="${TMP_MMORE_OUT}/mmore.results.stderr"
			MMORE_M8OUT="${TMP_MMORE_OUT}/mmore.results.m8out"
			MMORE_MYOUT="${TMP_MMORE_OUT}/mmore.results.myout"
			MMORE_MYDOMOUT="${TMP_MMORE_OUT}/mmore.results.mydomout"
			MMORE_MYTIME="${TMP_MMORE_OUT}/mmore.results.mytimeout"
			MMORE_MYTHRESH="${TMP_MMORE_OUT}/mmore.results.mythreshout"
		}

		# If overwriting is not permitted, check that destination files do not already exist
		{
			if [[ $DO_OVERWRITE == "0" ]];
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

		# Soft link input files into temporary working folders 
		{
			:
		}
	}

	# === GET DATABASE SIZES === #
	{
		# get target database size
		TARGET_DBSIZE=$( grep "^HMMER" ${TARGET_MMORE} | wc -l )
		echo_v 3 "TARGET_DBSIZE: $TARGET_DBSIZE"

		# get query database size
		QUERY_DBSIZE=$( grep "^>" ${QUERY_MMORE} | wc -l )
		echo_v 3 "QUERY_DBSIZE: $QUERY_DBSIZE"

		# get e-values from p-values 
		MMSEQS_EVAL=$( python <<< "eval = ${MMSEQS_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )
		MMSEQS_P2S_EVAL=$( python <<< "eval = ${MMSEQS_P2S_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )
		MMSEQS_S2S_EVAL=$( python <<< "eval = ${MMSEQS_S2S_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )
		MMORE_VITERBI_EVAL=$( python <<< "eval = ${MMORE_VITERBI_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )
		MMORE_CLOUD_EVAL=$( python <<< "eval = ${MMORE_CLOUD_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )
		MMORE_BOUNDFWD_EVAL=$( python <<< "eval = ${MMORE_BOUNDFWD_PVAL} * ${QUERY_DBSIZE}; print('%.3e' % eval)" )

		echo_v 3 "MMSEQS_EVAL: $MMSEQS_EVAL, MMSEQS_PVAL: $MMSEQS_PVAL"
	}

	# === MMSEQS SEARCH === #
	{
		echo "# RUN_MMSEQS?: $MMSEQS_DO_MMSEQS"

		if (( $MMSEQS_DO_MMSEQS != 0 ));
		then
		{
			# Configure MMSEQS Options
			{
				if (( $MMORE_DO_MMSEQS_ALN == 1 )); then
					MMSEQS_ALN="-a"
				else 
					MMSEQS_ALN=""
				fi 

				FORMAT_OUTPUT="qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
				# FORMAT_OUTPUT="qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar"
			}

			# Run MMSEQS prefilter
			{
				echo_v 3 "# Running: MMSEQS P2S Prefilter..." 
				BEG_TIME=$(GET_TIME)

				# determine whether to use sensitivity or kscore search params
				echo "MMSEQS_DO_SENS: $MMSEQS_DO_SENS"
				if (( $MMSEQS_DO_SENS == 1 ))
				then
					MMSEQS_KSCORE=""
				else
					MMSEQS_SENS=""
				fi 

				ECHO_AND_RUN 3 \
				$MMSEQS prefilter 															\
				$TARGET_MMSEQS_P $QUERY_MMSEQS											\
				$MMSEQS_PREFILTER 															\
																									\
				$(SET_ARGOPT 1 	-v 							$VERBOSE) 				\
				$(SET_ARGOPT 1 	--threads 					$NUM_THREADS) 			\
				$(SET_ARGOPT 1 	--max-seqs 					$MMSEQS_MAXSEQS) 		\
				$(SET_ARGOPT 1 	--k-score 					$MMSEQS_KSCORE) 		\
				$(SET_ARGOPT 1 	-s 							$MMSEQS_SENS) 			\
				$(SET_ARGOPT 1 	--min-ungapped-score 	$MMSEQS_UNGAPPED) 	\
				# > $MMSEQS_PREFILTER_STDOUT 												\

				CHECK_ERROR_CODE "PREFILTER"
				MMSEQS_PRV=$MMSEQS_PREFILTER

				END_TIME=$(GET_TIME)
				TIME_MMSEQS_PREFILTER=$(GET_DURATION $BEG_TIME $END_TIME)
				printf "# TIME_MMSEQS_P2S_PREFILTER: %.3f sec\n" $TIME_MMSEQS_PREFILTER
			}

			# If doing stats, find number of results
			if (( $DO_STATS == 1 ));
			then
			{
				PREFILTER_DBSIZE="0"
				DB_FILES=$( ls ${MMSEQS_PREFILTER}.[1-9]* )
				for FILE in $DB_FILES
				do
					DBSIZE=$( wc -l ${FILE} | awk '{ print $1 }' )
					PREFILTER_DBSIZE=$(( ${PREFILTER_DBSIZE} + ${DBSIZE} ))
				done	

				echo "# PREFILTER_DBSIZE: $PREFILTER_DBSIZE"
			}
			fi

			# If Profile-to-Sequence search
			if [[ $SEARCH_TYPE == "P2S" ]]
			then 
			{
				# Run MMSEQS Profile-to-Sequence alignment search
				{
					echo_v 3 "# Running: MMSEQS P2S Viterbi Search..."
					BEG_TIME=$(GET_TIME)

					ECHO_AND_RUN 3 \
					$MMSEQS align 															\
					$TARGET_MMSEQS_P  $QUERY_MMSEQS 									\
					$MMSEQS_PRV 															\
					$MMSEQS_P2S 															\
																								\
					$(SET_ARGOPT 1 	-v 				$VERBOSE) 					\
					$(SET_ARGOPT 1 	--threads 		$NUM_THREADS) 				\
					$(SET_ARGOPT 1 	-e					$MMSEQS_EVAL) 				\
					# > $MMSEQS_P2S_STDOUT												\

					CHECK_ERROR_CODE "P2S_ALIGN"
					MMSEQS_PRV=$MMSEQS_P2S

					END_TIME=$(GET_TIME)
					TIME_MMSEQS_P2S=$(GET_DURATION $BEG_TIME $END_TIME)
					printf "# TIME_MMSEQS_P2S_ALIGN: %.3f sec\n" $TIME_MMSEQS_P2S
				}

				# Report results to m8 file
				{
					echo_v 3 "# Running: MMSEQS P2S Convert Alignments..."
					BEG_TIME=$(GET_TIME)

					ECHO_AND_RUN 3 \
					$MMSEQS convertalis 										\
					$TARGET_MMSEQS_P  $QUERY_MMSEQS 						\
					$MMSEQS_P2S 												\
					$MMSEQS_P2S_M8 	 										\
																					\
					$(SET_ARGOPT 1 	-v 				$VERBOSE) 		\
					# > $MMSEQS_P2S_CONVERTALIS_REPORT					\

					CHECK_ERROR_CODE "CONVERT_ALIGN"

					END_TIME=$(GET_TIME)
					TIME_MMSEQS_CONVERTALIS=$(GET_DURATION $BEG_TIME $END_TIME)
					printf "# TIME_MMSEQS_P2S_CONVERTALIS: %.3f sec\n" $TIME_MMSEQS_CONVERTALIS
				}

				# If doing stats, find number of results
				{
					P2S_DBSIZE="0"
					DB_FILES=$( ls ${MMSEQS_P2S}.[1-9]* )
					for FILE in $DB_FILES
					do
						DBSIZE=$( wc -l ${FILE} | awk '{ print $1 }' )
						P2S_DBSIZE=$(( ${P2S_DBSIZE} + ${DBSIZE} ))
					done	

					echo "# P2S_DBSIZE: $P2S_DBSIZE"
				}

				# Translate MMSEQS Profile IDs to MMSEQS Sequence IDs
				{
					echo_v 3 "# Running: Translating P2S to S2S IDs..."
					BEG_TIME=$(GET_TIME)
					TRANSLATE_VERBOSE="0"

					ECHO_AND_RUN 3 \
					bash ${SCRIPT_DIR}/helpers/translate_ids.sh  	\
					$TARGET_MMSEQS_P 	$TARGET_MMSEQS_S 					\
					$MMSEQS_P2S 												\
					$MMSEQS_DUMMY 												\
					$TRANSLATE_VERBOSE 										\

					CHECK_ERROR_CODE "TRANSLATE_IDS"

					END_TIME=$(GET_TIME)
					TIME_MMSEQS_IDS=$(GET_DURATION $BEG_TIME $END_TIME)
					printf "# TIME_MMSEQS_TRANSLATE_IDS: %.3f sec\n" $TIME_MMSEQS_IDS
				}
			}
			fi

			# Run MMSEQS Sequence-to-Sequence search
			{
				echo_v 3 "# Running: MMSEQS S2S Search..."
				BEG_TIME=$(GET_TIME)

				ECHO_AND_RUN 3 \
				$MMSEQS align 																\
				$TARGET_MMSEQS_S  $QUERY_MMSEQS 										\
				$MMSEQS_PRV 																\
				$MMSEQS_S2S 																\
																								\
				$(SET_ARGOPT 1 	-v 					$VERBOSE) 					\
				$(SET_ARGOPT 1 	--threads 			$NUM_THREADS) 				\
				$(SET_ARGOPT 1 	-e 					"inf" )						\
				$(SET_ARGOPT 1 	--alt-ali			"0" ) 						\
				$(SET_ARGOPT 1 	-a 					"0" ) 						\
				$(SET_ARGOPT 1 	--alignment-mode	"3" ) 						\
				$(SET_ARGOPT 1 	--realign 			"0" ) 						\
				# > $MMSEQS_S2S_REPORT 												 	\

				CHECK_ERROR_CODE "S2S_ALIGN"

				END_TIME=$(GET_TIME)
				TIME_MMSEQS_S2S=$(GET_DURATION $BEG_TIME $END_TIME)
				printf "# TIME_MMSEQS_S2S_ALIGN: %.3f sec\n" $TIME_MMSEQS_S2S
			}

			# If doing stats, find number of results
			{
				S2S_DBSIZE="0"
				DB_FILES=$( ls ${MMSEQS_S2S}.[1-9]* )
				for FILE in $DB_FILES
				do
					DBSIZE=$( wc -l ${FILE} | awk '{ print $1 }' )
					S2S_DBSIZE=$(( ${S2S_DBSIZE} + ${DBSIZE} ))
				done	

				printf "# S2S_DBSIZE: $S2S_DBSIZE"
			}

			# Report results to m8 file
			{
				echo_v 3 "# Running: MMSEQS S2S Convert Alignments..."
				BEG_TIME=$(GET_TIME)

				ECHO_AND_RUN 3 \
				$MMSEQS convertalis 										\
				$TARGET_MMSEQS_S  $QUERY_MMSEQS 						\
				$MMSEQS_S2S 												\
				$MMSEQS_S2S_M8 	 										\
																				\
				$(SET_ARGOPT 1 	-v 				$VERBOSE) 		\
				$(SET_ARGOPT 1 	--threads 		$NUM_THREADS) 	\
				# > $MMSEQS_S2S_CONVERTALIS_REPORT					\

				CHECK_ERROR_CODE "CONVERT_ALIGN"

				END_TIME=$(GET_TIME)
				TIME_MMSEQS_CONVERTALIS=$(GET_DURATION $BEG_TIME $END_TIME)
				printf "# TIME_MMSEQS_S2S_CONVERTALIS: %.3f sec\n" $TIME_MMSEQS_CONVERTALIS
			}

			# Add P2S Viterbi E-values from P2S Search
			if [[ $SEARCH_TYPE == "P2S" ]]
			then 
			{
				echo_v 3 "# Sorting P2S Results..."
				# sort P2S results
				cat ${MMSEQS_P2S_M8} | sort -k1,1 -k2,2 -s > ${MMSEQS_P2S_M8}.sort
				mv ${MMSEQS_P2S_M8}.sort ${MMSEQS_P2S_M8}
				# extract P2S fields
				# cat ${MMSEQS_P2S_M8} | awk '{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' > ${MMSEQS_P2S_M8}.align
				cat ${MMSEQS_P2S_M8} | cut -f 11-12 > ${MMSEQS_P2S_M8}.eval
				
				echo_v 3 "# Sorting S2S Results..."
				# sort S2S results
				cat ${MMSEQS_S2S_M8} | sort -k1,1 -k2,2 -s > ${MMSEQS_S2S_M8}.sort
				mv ${MMSEQS_S2S_M8}.sort ${MMSEQS_S2S_M8}
				# extract S2S fields
				cat ${MMSEQS_S2S_M8} | cut -f 1-10 > ${MMSEQS_S2S_M8}.align
				# cat ${MMSEQS_S2S_M8} | awk '{ print $11,$12 }' > ${MMSEQS_S2S_M8}.eval
				
				echo_v 3 "# Merging S2S and P2S Results..."
				# merge results
				paste ${MMSEQS_S2S_M8}.align ${MMSEQS_P2S_M8}.eval > ${MMSEQS_M8}

				# remove temporary files
				# rm ${MMSEQS_P2S_M8}.align 
				rm ${MMSEQS_P2S_M8}.eval
				rm ${MMSEQS_S2S_M8}.align 
				# rm ${MMSEQS_S2S_M8}.eval
			}
			else 
			{
				:
			}
			fi 

			# Copy results out to output location
			IF_THEN_COPY "${MMSEQS_M8}" "${MMSEQS_M8_OUT}" "${MMSEQS_M8_OUT}"
		}
		else 
		{
			echo "# WARNING: MMSEQS flagged not to run [ --run-mmseqs 0 ]"
		}
		fi

	}

	# === CONVERT MMSEQS OUTPUT -> MMORE INPUT === #
	{
		:
	}

	# === MMORE SEARCH === #
	{
		TARGET_INDEX="${TMP_MMORE}/target.idx" 
		QUERY_INDEX="${TMP_MMORE}/query.idx"

		echo "# RUN_MMORE?: $MMORE_DO_MMORE"

		# Run Search
		if (( "MMORE_DO_MMORE" != "0" ))
		then
		{

			ECHO_AND_RUN 3 \
			$MMORESEQS  mmore-search 													\
			$TARGET_MMORE 	$QUERY_MMORE 												\
			$MMSEQS_M8	 																	\
																								\
			$(SET_ARGOPT 1 	--program-mmoreseqs 	$MMORESEQS) 				\
			$(SET_ARGOPT 1 	--program-mmseqs 		$MMSEQS) 					\
			$(SET_ARGOPT 1 	--program-hmmer 		$HMMBUILD) 					\
																								\
			$(SET_ARGOPT 1 	--verbose 				$VERBOSE) 					\
			$(SET_ARGOPT 1 	--num-threads 			$NUM_THREADS) 				\
																								\
			$(SET_ARGOPT 1 	--run-full 				$MMORE_DO_FULL) 			\
			$(SET_ARGOPT 1 	--run-bias 				$MMORE_DO_DOMAIN) 		\
			$(SET_ARGOPT 1 	--run-mmseqsaln 		$MMORE_DO_MMSEQS_ALN) 	\
			$(SET_ARGOPT 1 	--run-vitaln 			$MMORE_DO_VIT_ALN) 		\
			$(SET_ARGOPT 1    --run-postaln 			$MMORE_DO_POST_ALN) 		\
																								\
			$(SET_ARGOPT 2 	--dbsizes 				$TARGET_DBSIZE   			\
																$QUERY_DBSIZE) 			\
																								\
			$(SET_ARGOPT 1 	--alpha 					$MMORE_ALPHA) 				\
			$(SET_ARGOPT 1 	--beta 					$MMORE_BETA) 				\
			$(SET_ARGOPT 1 	--gamma 					$MMORE_GAMMA) 				\
																								\
			$(SET_ARGOPT 1 	--run-filter 			$MMORE_DO_FILTER) 		\
			$(SET_ARGOPT 1 	--run-vit-filter 		$MMORE_DO_VIT_FILTER) 	\
			$(SET_ARGOPT 1 	--run-cld-filter 		$MMORE_DO_CLD_FILTER) 	\
			$(SET_ARGOPT 1 	--run-fwd-filter 		$MMORE_DO_FWD_FILTER) 	\
			$(SET_ARGOPT 1 	--filter-vit 			$MMORE_VITERBI) 			\
			$(SET_ARGOPT 1 	--filter-cld 			$MMORE_CLOUD) 				\
			$(SET_ARGOPT 1 	--filter-fwd 			$MMORE_FWDBACK) 			\
																								\
			$(SET_ARGOPT 1 	--m8out 					$MMORE_M8OUT)				\
			$(SET_ARGOPT 1 	--myout 					$MMORE_MYOUT)				\
			$(SET_ARGOPT 1 	--mydomtblout 			$MMORE_MYDOMOUT)			\
			$(SET_ARGOPT 1 	--mytimeout 			$MMORE_MYTIME)				\
			$(SET_ARGOPT 1 	--mythreshout 			$MMORE_MYTHRESH)			\

			CHECK_ERROR_CODE "MMORE_MAIN"
		}
		else 
		{
			echo "# WARNING: MMORE flagged not to run [ --run-mmore 0 ]"
		}
		fi

		# Copy results out to output location

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

	local END_TOTAL=$(GET_TIME)
	TIME_PREP=$(GET_DURATION $BEG_TOTAL $END_TOTAL)
	printf "# SEARCH_TIME: %.3f sec\n" $TIME_PREP
}
