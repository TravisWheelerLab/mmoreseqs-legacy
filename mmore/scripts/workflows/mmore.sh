#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore.sh	
#	AUTHOR:		David Rich
#	DESC: 		runs mmseqs, then pipes output to mmore search
###########################################################################

# PROGRAM:
# (1) mmseqs creates a temporary copy of query/target dbs
# (2) mmseqs creates index of query/target dbs
# (3) mmore creates index of query/target dbs
# (4) mmseqs search
# (5) reformat mmseqs results for input to mmore
# (5) mmore search

# only echoes if has higher verbosity level
function echo_v
{
	if (( $VERBOSE >= $1 ))
	then
		echo $2
	fi
}
VERBOSE="${VERBOSE:-3}"
echo "VERBOSE LEVEL: $VERBOSE"

# commandline variables (verify proper number of variables)
NUM_ARGS=$#
if (( $NUM_ARGS < 3 )); then 
	echo "ERROR: illegal number of parameters: $NUM_ARGS of 3"
	echo "Usage: <target> <query> <target_mmseqs>"
	exit
fi

ARG_TARGET=$1
ARG_QUERY=$2
ARG_TARGET_MMSEQS=$3

echo_v 1 "MMORE SCRIPT: $ARG_TARGET $ARG_QUERY $ARG_TARGET_MMSEQS"

# if system tools allowed, then try to using system installs
USE_LOCAL_TOOLS="${USE_LOCAL_TOOLS:-0}"
if (( USE_LOCAL_TOOLS == 0 ))
then
	echo_v 3 "SEARCHING FOR SYSTEM TOOLS..."
	# which mmseqs
	MMSEQS_DEFAULT=$(which mmseqs)
	# which hmmbuild
	HMMER_DEFAULT=$(which hmmbuild)
	# which mmore
	MMORE_DEFAULT=$(which mmore)
fi

# if not install on system, add local install to path (paths passed as environmental vars in main)
if [ -z "$MMSEQS_DEFAULT" ]
then
	echo_v 3 "USING LOCAL MMSEQS: ${MMSEQS_DIR}"
	export PATH="${MMSEQS_DIR}:$PATH"
fi
if [ -z "$HMMER_DEFAULT" ]
then
	echo_v 3 "USING LOCAL HMMER: ${HMMER_DIR}"
	export PATH="${HMMER_DIR}:$PATH"
fi
if [ -z "$MMORE_DEFAULT" ]
then
	echo_v 3 "USING LOCAL MMORE: ${MMORE_DIR}"
	export PATH="${MMORE_DIR}:$PATH"
fi

# if program still can't be found, throw error.
# which mmseqs
MMSEQS=$(which mmseqs)
if [ -z "$MMSEQS" ]
then
	echo "ERROR: hmmer is not installed on system or locally."
	exit 1
fi
# which hmmbuild
HMMBUILD=$(which hmmbuild)
if [ -z "$HMMBUILD" ]
then
	echo "ERROR: hmmer is not installed on system or locally."
	exit 1
fi
# which mmore
MMORE=$(which mmore)
if [ -z "$MMORE" ]
then
	echo "ERROR: mmore is not installed on system or locally."
	exit 1
fi

# if args not supplied, falls back on these defaults
# main args:
TARGET="${TARGET:-$TARGET_ARG}"
QUERY="${QUERY:-$QUERY_ARG}"
TARGET_MMSEQS="${TARGET_MMSEQS:-$TARGET_MMSEQS_ARG}"
# main arg types:
TARGET_TYPE="${TARGET_TYPE:-'HMM'}"
QUERY_TYPE="${QUERY_TYPE:-'FASTA'}"
TARGET_MMSEQS_TYPE="${TARGET_MMSEQS_TYPE:-'HHM'}"
# options :
ROOT_DIR="${ROOT_DIR:-'/'}"
TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-XXXX)}"
RM_TEMP="${RM_TEMP:-0}" 
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

# output files:

# report variables
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

# ##### MAIN ###### #

BENCH_DIR=$(pwd)

# temporary directory
TMP_ROOT=${TEMP_DIR}/
TMP=${TMP_ROOT}/

# temporaty subdirectories for mmseqs and 
TMP_MMSEQS=${TMP}/mmseqs/
TMP_CLOUD=${TMP}/cloud/

# temporary databases for mmseqs
QUERY_MMSEQS=${TMP_MMSEQS}/query
TARGET_MMSEQS=${TMP_MMSEQS}/target

# result destination files
RESULT_RAW_MMSEQS=${TMP_MMSEQS}/results
RESULT_MMSEQS=${TMP_MMSEQS}/results.m8
RESULT_PLUS_MMSEQS=${TMP_MMSEQS}/results.m8+
RESULT_CLOUD=${TMP_CLOUD}/results.m8

# make tmp directories
mkdir $TMP_ROOT
mkdir $TMP
mkdir $TMP_CLOUD
mkdir $TMP_MMSEQS

# mmseqs options
MMSEQS_SHUFFLE=""

# === FILE FORMATTING === #
echo "# (0/4) File formatting and building..."

# If query is a hmm file, then we need to 
# build consensus fasta file for all hmms in database
# WARNING: this may not be intended
if [[ "$QUERY_TYPE" = "HMM" ]]
then
	QUERY_HMM=$QUERY
	echo "WARNING: Query is an HMM file. Profile-to-Profile is not supported, \
			so query being converted to FASTA consensus sequences. 				  \
			This will result in a loss of accuracy."
fi

# if target mmseqs file is a hhm, then we are ready.
if [[ "$QUERY_TYPE" = "FASTA" ]]
then
	QUERY_FASTA=$QUERY
fi

# Final check if file is a valid file format
if [[ "$QUERY_TYPE" != "FASTA" && "$QUERY_TYPE" != "HMM" ]]
then
	echo "ERROR: Query is an unsupported filetype of '$QUERY_TYPE'."
	exit 1
fi

# if target is a msa file, then we need to
# generate an HMM file from it using hmmer suite
if [[ "TARGET_TYPE" = "MSA" ]]
then 
	TARGET_MSA=$TARGET 

fi

# If target is a hmm file, then we need to 
# build consensus fasta file for all hmms in database
if [[ "$TARGET_TYPE" != "HMM" ]]
then
	echo "ERROR: Target is an unsupported filetype of '$TARGET_TYPE'."
	exit 1
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

# if target mmseqs file is a msa file, then we need to
# convert it to mm_msa, then generate 
if [[ "$TARGET_MMSEQS_TYPE" == "MSA" ]]
	TARGET_MMSEQS_MSA=$TARGET_MMSEQS
fi

# if target mmseqs file is a msa file, then we need to
# convert it to mm_msa, then generate 
if [[ "$TARGET_MMSEQS_TYPE" == "MM_MSA" ]]
	TARGET_MMSEQS_MM_MSA=$TARGET_MMSEQS
fi

# if target mmseqs file is a hhm, then we are ready.
if [[ "$TARGET_MMSEQS_TYPE" == "HHM" ]]
	TARGET_MMSEQS_HHM=$TARGET_MMSEQS
fi

# Final check if file is a valid file format
if [[ "$TARGET_MMSEQS_TYPE" != "HHM" && "$TARGET_MMSEQS_TYPE" != "FASTA" ]]
	echo "ERROR: MMseqs target is an unsupported type of '$TARGET_MMSEQS_TYPE'."
	exit 1
fi

exit

# ==== CLOUD INDEX ==== #
echo "# (1/4) Create CLOUD index of database..."
time $CLOUD index $TARGET $QUERY
time mv $TARGET.idx $TMP_CLOUD/target.idx 
time mv $QUERY.idx $TMP_CLOUD/query.idx

# capture database size
NUM_TARGETS=$( grep "#" -v $TMP_CLOUD/target.idx | wc -l )
NUM_QUERIES=$( grep "#" -v $TMP_CLOUD/query.idx | wc -l )
DB_SIZE=$((NUM_TARGETS))

# ==== MMSEQS SEARCH ======= $
echo "# (2/4) Performing MMSEQS Search..."
# --format-output STR       	CSV list of output columns from: query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid,taxid,taxname,taxlineage  [query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits]
# -s FLOAT                     	Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [5.700]
# -k INT                       	k-mer length (0: automatically set to optimum) [0]
# --k-score INT                	k-mer threshold for generating similar k-mer lists [2147483647]
# --max-rejected INT           	Maximum rejected alignments before alignment calculation for a query is stopped [2147483647]
# --max-accept INT             	Maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
# --local-tmp STR               Path where some of the temporary files will be created []
# --remove-tmp-files BOOL       Delete temporary files [1]
# --e-profile FLOAT             Include sequences matches with < e-value thr. into the profile (>=0.0) [0.001]

FORMAT_OUTPUT="qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
EVAL_CUTOFF=$(( PVAL_CUTOFF * DB_SIZE ))
EVAL_REPORT=$(( PVAL_REPORT * DB_SIZE ))

time $MMSEQS easy-search 
	$QUERY_MMSEQS  	$TARGET_MMSEQS 						\
	$RESULT_MMSEQS  	$TMP_MMSEQS 							\
	-k 							$KMER 							\
	--min-ungapped-score 	$MIN_UNGAPPED_SCORE 			\
	-e 							$EVAL_CUTOFF 					\
	--k-score 					$K_SCORE 						\
	--remove-tmp-files  		0									\

# ==== FBPRUNER SEARCH ==== #
echo "(7/7) Performing FB_PRUNER Search..."
# cloud search options:
# --alpha FLOAT 				alpha cloud tuning parameter
# --beta INT 					beta cloud tuning parameter
# --pipeline STR/INT 			Choose: [test,main,mmseqs,index]
# --threshold FLOAT 			Score Reporting Threshold
# --format-input STR 			CSV list of input format
# --mmseqs-results STR 			Location of mmseqs results
# --mmseqs-tid STR 				Location of mmseqs target id lookup table
# --mmseqs-qid STR 				Location of mmseqs query id lookup table 

TARGET_INDEX="$TMP_CLOUD/target.idx" 
QUERY_INDEX="$TMP_CLOUD/query.idx"

# check if output are set
if [ ! -z "$STDOUT" ]; then
	STDOUT=""
else
	STDOUT="--stdout $STDOUT"
fi
if [ ! -z "$M8_OUT" ]; then
	M8_OUT=""
else
	$STDOUT="--m8-out $STDOUT"
fi
if [ ! -z "$STDOUT" ]; then
	$STDOUT=""
else
	$STDOUT="--stdout $STDOUT"
fi


time $FB_PRUNER mmore_main 								\
	$TARGET_PATH $QUERY_PATH 								\
	--mmseqs-results 	$RESULT_MMSEQS						\
	--index 			$TARGET_INDEX $QUERY_INDEX 		\
	$STDOUT 														\
	$TBL_OUTPUT													\
	$M8_OUTPUT 													\
	$MY_OUTPUT													\
	--alpha 				$ALPHA 								\
	--beta 				$BETA 								\
	--gamma 				$GAMMA 								\

mv $RESULT_CLOUD $RESULTS

# === CLEAN UP === #
# remove temporary folders
if ( $REMOVE_TEMP > 0 )
then
	rm -r $TEMP_DIR
fi
