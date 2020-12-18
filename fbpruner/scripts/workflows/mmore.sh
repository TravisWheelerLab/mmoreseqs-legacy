#!/usr/bin/bash
###########################################################################
#	NAME: 		mmseqs_plus.sh	
#	AUTHOR:		David Rich
#	DESC: 		runs mmseqs, then pipes output to fb-pruner search
###########################################################################

# PROGRAM:
# (1) mmseqs creates a temporary copy of query/target dbs
# (2) mmseqs creates index of query/target dbs
# (3) fb-pruner creates index of query/target dbs
# (4) mmseqs search
# (5) reformat mmseqs results for input in fb-pruner
# (5) fb-pruner search

# SET ENVIRONMENTAL VARIABLES

# commandline variables (verify proper number of variables)
NUM_ARGS=$#
if (( $NUM_ARGS < 25 )); then 
	echo "ERROR: illegal number of parameters: $NUM_ARGS of 25"
	echo "Usage: <target> <query> <results>"
	exit
fi

# -- Main Args -- #
TARGET=$1
QUERY=$2
T_TYPE=$3
Q_TYPE=$4
# -- Options -- #
RESULTS=$5
RM_TEMP=$6
VERBOSE=$7
# -- MMseqs Args -- #
KMER=$8
KSCORE=$9
MIN_UNGAPPED_SCORE=$10
PVAL_CUTOFF=$11
EVAL_CUTOFF=$12
# -- MMORE-SEQS / FB-PRUNER Args -- #
ALPHA=$13
BETA=$14
GAMMA=$15
PVAL_REPORT=$16
EVAL_REPORT=$17
COMP_BIO=$18
# -- TOOL LOCATION -- #
FASTA_T0_HMM_SCRIPT=$19
MMSEQS_BIN=$20
HMMBUILD_BIN=$21
# -- OUTPUTS -- #
STDOUT=$22
TBL_OUT=$23
M8_OUT=$24
MY_OUT=$25

# default parameters
ROOT_DIR="${ROOT_DIR:-'/'}"
TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-XXXX)}"
RM_TEMP="${RM_TEMP:-0}" 
# mmseqs parameters
KMER="${KMER:-7}"
K_SCORE="${K_SCORE:-75}"
MIN_UNGAPPED_SCORE="${MIN_UNGAPPED_SCORE:-15}"
PVAL_CUTOFF="${PVAL_CUTOFF:-0.001}"
EVAL_CUTOFF="${EVAL_CUTOFF:-200}"
SENSITIVITY="${SENSITIVITY:-7.5}"
# fbpruner
ALPHA="${ALPHA:-12.0}"
BETA="${BETA:-16.0}"
GAMMA="${GAMMA:-5}"					
BIAS_FILTER="${BIAS_FILTER:-1}"  	# boolean whether to compute bias filter
EVAL_REPORT="${EVAL_REPORT:-0.005}"   
PVAL_REPORT="${PVAL_REPORT:-}"

# report variables
if ( $VERBOSE > 3 ); then
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

exit

# programs
MMSEQS=mmseqs
CLOUD=./build/fb-pruner

# main
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
	$QUERY_MMSEQS  $TARGET_MMSEQS 					\
	$RESULT_MMSEQS  $TMP_MMSEQS 					\
	-k 						$KMER 					\
	--min-ungapped-score 	$MIN_UNGAPPED_SCORE 	\
	-e 						$EVAL_CUTOFF 			\
	--k-score 				$K_SCORE 				\
	--remove-tmp-files  	0						\

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
if [ "$STDOUT" = "==" ]; then
	STDOUT=""
else
	STDOUT="--stdout $STDOUT"
fi
if [ "$M8_OUT" = "==" ]; then
	M8_OUT=""
else
	$STDOUT="--m8-out $STDOUT"
fi
if [ "$STDOUT" = "==" ]; then
	$STDOUT=""
else
	$STDOUT="--stdout $STDOUT"
fi


time $FB_PRUNER mmseqs 								\
	$TARGET_PATH $QUERY_PATH 						\
	--mmseqs-results 	$RESULT_MMSEQS				\
	--index 			$TARGET_INDEX $QUERY_INDEX 	\
	$STDOUT 										\
	$TBL_OUTPUT										\
	$M8_OUTPUT 										\
	$MY_OUTPUT										\
	--alpha 			$ALPHA 						\
	--beta 				$BETA 						\
	--gamma 			$GAMMA 						\

mv $RESULT_CLOUD $RESULTS

# === CLEAN UP === #
# remove temporary folders
if ( $REMOVE_TEMP > 0 )
then
	rm -r $TEMP_DIR
fi
