#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore-prep.sh	
#	AUTHOR:		David Rich
#	DESC: 		Creates any necessary files to prepare to run main pipeline.
###########################################################################

# PROGRAM:
# (1a) if <target> is .msa file:
# (1b) convert <target> to .hmm file
# (1b) convert <target> to .hmm file
#
# (2a) if <query> is .msa file, convert to .hmm file
# (2b) if <query> is .hmm file, convert to .fasta (consensus) file
# (2c) verify that final <query> is .fasta file
#

# ##### FUNCTIONS ###### #
# main script only contains the minimum functions needed to load helper function script
{
	#set verbosity (MAXIMUM if not specified)
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

	# load script dependencies 
	{
		# load functions
		echo_v 3 "Importing 'mmore_functions.sh..."
		LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_functions.sh"
		# load tools 
		echo_v 3 "Importing 'mmore_get-tools.sh'..."
		LOAD_SOURCE "${SCRIPT_DIR}/helpers/mmore_get-tools.sh"
	}
}

# ##### MAIN ###### #
{
	echo_v 2 "# Running Program: mmore-prep.sh..."
	PROGRAM+="> mmore-prep.sh >"

	BEG_TOTAL=$(GET_TIME)

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
			ARG_PREP_DIR=$3
			ARG_TARGET_TYPE=$4
			ARG_QUERY_TYPE=$5
		}

		# process environmental variables: if optional args not supplied by environment, falls back on these defaults
		{
			# prep args:
			TARGET="${TARGET:-$ARG_TARGET}"
			QUERY="${QUERY:-$ARG_QUERY}"
			TEMP_DIR="${TEMP_DIR:-$ARG_TEMP_DIR}"
			PREP_DIR="${PREP_DIR:-$ARG_PREP_DIR}"
			TARGET_TYPE="${TARGET_TYPE:-$ARG_TARGET_TYPE}"
			QUERY_TYPE="${QUERY_TYPE:-$ARG_QUERY_TYPE}"
			
			# load defaults:
			SET_ENV_ARG_DEFAULTS
		}

		# report variables
		{
			if (( $VERBOSE >= 3 )); then
				echo "# ========== MMORE: PREP ===================="
				echo "#      TARGET_PREP:  $TARGET [$TARGET_TYPE]"
				echo "#       QUERY_PREP:  $QUERY [$QUERY_TYPE]"
				echo "#         PREP_DIR:  $TEMP_DIR"
				echo "# ==========================================="
			fi
		}

		# build temporary folders
		{
			# top-level temporary directory
			TMP="${TEMP_DIR}/"
			TMP_DB="${TMP}/db/"
			TMP_OUT="${TMP}/out/"
			# mmore folders
			TMP_MMORE="${TMP}/mmore/"
			TMP_MMORE_DB="${TMP_MMORE}/db/"
			TMP_MMORE_OUT="${TMP_MMORE}/out/"
			# mmseqs folders
			TMP_MMSEQS="${TMP}/mmseqs/"
			TMP_MMSEQS_DB="${TMP_MMSEQS}/db/"
			TMP_MMSEQS_OUT="${TMP_MMSEQS}/out/"

			# make tmp directories
			MAKE_DIR "$TMP"
			MAKE_DIR "$TMP_DB"
			MAKE_DIR "$TMP_OUT"
			MAKE_DIR "$TMP_MMORE"
			MAKE_DIR "$TMP_MMORE_DB"
			MAKE_DIR "$TMP_MMORE_OUT"
			MAKE_DIR "$TMP_MMSEQS"
			MAKE_DIR "$TMP_MMSEQS_DB"
			MAKE_DIR "$TMP_MMSEQS_OUT"

			# list of temp folders 
			TEMP_FOLDERS=""
			TEMP_FOLDERS+="${TMP}" 
			TEMP_FOLDERS+="${TMP_DB}"
			TEMP_FOLDERS+="${TMP_OUT}"
			TEMP_FOLDERS+="${TMP_MMORE}"
			TEMP_FOLDERS+="${TMP_MMORE_DB}"
			TEMP_FOLDERS+="${TMP_MMORE_OUT}"
			TEMP_FOLDERS+="${TMP_MMSEQS}" 
			TEMP_FOLDERS+="${TMP_MMSEQS_DB}"
			TEMP_FOLDERS+="${TMP_MMSEQS_OUT}"
		}
	}

	# file formatting
	{
		# set whether to out
		if (( $VERBOSE < 3 )); then
			QUIET_OUT="/dev/null"
		fi
		
		# destination locations 
		{
			# main db
			TARGET_MSA="${TMP_DB}/target.msa"
			TARGET_HMM="${TMP_DB}/target.hmm"
			TARGET_FASTA="${TMP_DB}/target.fasta"
			# mmore db (soft link to main db)
			TARGET_MMORE="${TMP_MMORE_DB}/target.hmm"
			# mmseqs db
			TARGET_MM_MSA="${TMP_MMSEQS_DB}/target.mm_msa"
			TARGET_MMSEQS_S="${TMP_MMSEQS_DB}/target.smmdb"
			TARGET_MMSEQS_P="${TMP_MMSEQS_DB}/target.pmmdb"

			# main db
			QUERY_MSA="${TMP_DB}/query.msa"
			QUERY_HMM="${TMP_DB}/query.hmm"
			QUERY_FASTA="${TMP_DB}/query.fasta"
			# mmore db (soft link to main db)
			QUERY_MMORE="${TMP_MMORE_DB}/query.fasta" 
			# mmseqs db
			QUERY_MMSEQS="${TMP_MMSEQS_DB}/query.smmdb"
		}
		
		# IMPORT main input files (MSA or FASTA files)
		{
			# import input files using copy or symbolic links
			{
				LINK="ln -nfs -v"
				COPY="cp -v"

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
				echo "ERROR: TARGET_TYPE '${TARGET_TYPE}' is not a valid file type."
				exit 1
			fi 

			# verify files exist 
			CHECK=$(CHECK_FILE_EXISTS ${TARGET})
			echo "CHECK: $TARGET $CHECK." 
			CHECK=$(CHECK_FILE_EXISTS ${QUERY})
			echo "CHECK: $QUERY $CHECK." 
			
			# import target file 
			echo "IMPORT: $IMPORT ${TARGET} $DO_COPY"
			TARGET_IN_TYPE="$TARGET_TYPE"
			if [ "$TARGET_TYPE" == "MSA" ]; then
				$IMPORT "${TARGET}" "${TARGET_MSA}"
			elif [ "$TARGET_TYPE" == "FASTA" ]; then
				$IMPORT "${TARGET}" "${TARGET_FASTA}" 
			fi 

			# verify valid query formats 
			if [ "$QUERY_TYPE" != "MSA" ] && [ "$QUERY_TYPE" != "FASTA" ]
			then 
				echo "ERROR: QUERY_TYPE '${QUERY_TYPE}' is not a valid file type."
				exit 1
			fi 

			# import query file
			QUERY_IN_TYPE="$QUERY_TYPE"
			if [ "$QUERY_TYPE" == "MSA" ]; then
				$IMPORT "${QUERY}" "${QUERY_MSA}"
			elif [ "$QUERY_TYPE" == "FASTA" ]; then
				$IMPORT "${QUERY}" "${QUERY_FASTA}"
			fi
		}

		echo_v 2 "# (1/3) Importing files: ${TARGET} [${TARGET_TYPE}], ${QUERY} [${QUERY_TYPE}]..."
		# === IMPORT FILES === #
		{
			:
		}

		echo_v 2 "# (2/3) Formatting TARGET [${TARGET_TYPE}] files..."
		# TARGET file formatting
		{
			TARGET_TYPE="$TARGET_IN_TYPE"

			# If target is a msa file 
			if [ "$TARGET_TYPE" == "MSA" ]
			then
				# search type
				SEARCH_TYPE="P2S"

				# convert msa to hmm file 
				echo_v 3 "# TARGET: MSA => HMM"
				$HMMBUILD "$TARGET_HMM" "$TARGET_MSA" 		\

				# convert hmm to fasta
				echo_v 3 "# TARGET: HMM => FASTA"
				$HMMEMIT -c "$TARGET_HMM" > "$TARGET_FASTA"	
				sed -i 's/-consensus//g' "$TARGET_FASTA"

				# link hmm to mmore
				echo_v 3 "# TARGET: HMM => MMORE"
				# $LINK $TARGET_HMM $TARGET_MMORE
				$COPY "$TARGET_HMM" "$TARGET_MMORE"

				# convert msa to mm_msa
				echo_v 3 "# TARGET: MSA => MM_MSA"
				$MMSEQS convertmsa								\
				"$TARGET_MSA" "$TARGET_MM_MSA" 				\
				--identifier-field	0							\
				-v 						$VERBOSE					\

				# convert fasta to sequence mmdb 
				echo_v 3 "# TARGET: FASTA => MMSEQS_S"
				$MMSEQS createdb 									\
				"$TARGET_FASTA" 	"$TARGET_MMSEQS_S" 		\
				-v 						$VERBOSE 				\

				# convert mm_msa to profile mmdb
				echo_v 3 "# TARGET: MM_MSA => MMSEQS_P"
				$MMSEQS msa2profile  							\
				"$TARGET_MM_MSA" "$TARGET_MMSEQS_P" 		\
				-v 						$VERBOSE 				\

			fi

			# If target is a hmm file, then convert to consensus fasta file
			if [ "$TARGET_TYPE" == "FASTA" ]
			then
				# search type
				SEARCH_TYPE="S2S"

				# convert fasta to hmm
				echo_v 3 "# TARGET: FASTA => HMM"
				${SCRIPT_DIR}/helpers/convert_fasta-to-hmm.sh 	\
				"$TARGET_FASTA" "$TARGET_HMM" "tmp_split" 		\

				# convert fasta to sequence mmdb 
				echo_v 3
				$MMSEQS createdb  										\
				"$TARGET_FASTA" "$TARGET_MMSEQS_S"  				\
				-v 						$VERBOSE 						\

			fi
		}

		echo_v 2 "# (3/3) Formatting QUERY [${QUERY_TYPE}] files..."
		# QUERY file formatting
		{
			# If query is a msa file
			if [ "$QUERY_IN_TYPE" == "MSA" ]
			then
				# convert msa to hmm
				$HMMBUILD  										\
				"$QUERY_HMM" "$QUERY_MSA" 					\

				# convert hmm to fasta
				$HMMEMIT -c "$QUERY_HMM" > "$QUERY_FASTA"	
				sed -i 's/-consensus//g' "$QUERY_FASTA"

				# link fasta to mmore
				# $LINK "$QUERY_FASTA" "$QUERY_MMORE"
				$COPY "$QUERY_FASTA" "$QUERY_MMORE"

				# convert fasta to sequence mmdb 
				$MMSEQS createdb  							\
				"$QUERY_FASTA" "$QUERY_MMSEQS"  			\

			fi

			# if query is a fasta file
			if [ "$QUERY_IN_TYPE" == "FASTA" ]
			then
				# link fasta to mmore
				# $LINK "$QUERY_FASTA" "$QUERY_MMORE"
				$COPY "$QUERY_FASTA" "$QUERY_MMORE"

				# convert fasta to sequence mmdb
				$MMSEQS createdb  							\
				"$QUERY_FASTA" "$QUERY_MMSEQS" 			\

			fi
		}

		echo_v 2 "# (4/4) Cleanup"
		# CLEANUP
		{
			SEARCH_FILE="${TMP}/search_type.txt"
			echo "$SEARCH_TYPE" > "$SEARCH_FILE"
		}
	}

	END_TOTAL=$(GET_TIME)
	TIME_PREP=$(GET_DURATION $BEG_TOTAL $END_TOTAL)
	printf "# TIME_PREP: %.3f sec\n" $TIME_PREP
}
