#!/usr/local/bin/bash
###########################################################################
#	NAME: 		mmseqs_plus.sh	
#	AUTHOR:		David Rich
#	DESC: 		runs mmseqs, then pipes output to fb-pruner search
###########################################################################

# NOTES: to reinstall mmseqs and update scripts => 
# $ cd build/ 
# $ make clean 
# $ make -j 4 
# $ make install

# PROGRAM:
# (1) mmseqs creates a temporary copy of query/target dbs
# (2) mmseqs creates index of query/target dbs
# (3) fb-pruner creates index of query/target dbs
# (4) mmseqs search
# (5) reformat mmseqs results for input in fb-pruner
# (5) fb-pruner search

# commandline variables (verify proper number of variables)
NUM_ARGS=$#
if (( $NUM_ARGS != 3 )); then 
	echo "ERROR: illegal number of parameters"
	echo "Usage: <target> <query> <results>"
	exit
fi
TARGET=$1
QUERY=$2
RESULTS=$3

# default parameters
ROOT_DIR="${ROOT_DIR:-'/'}"
TEMP_DIR="${TEMP_DIR:-$(mktemp -d tmp-XXXX)}"
RM_TEMP="${RM_TEMP:-0}" 
# mmseqs
PVAL_CUTOFF="${PVAL_CUTOFF:-0.001}"
K_SCORE="${K_SCORE:-75}"
KMER="${KMER:-7}"
SENSITIVITY="${SENSITIVITY:-7.5}"
# fbpruner
ALPHA="${ALPHA:-12.0}"
BETA="${BETA:-16.0}"
GAMMA="${GAMMA:-5}"


# report variables
echo "#  TARGET: $TARGET"
echo "#   QUERY: $QUERY"
echo "# RESULTS: $RESULTS"
echo "#   ALPHA: $ALPHA"
echo "#    BETA: $BETA"
echo "#   GAMMA: $GAMMA"

# programs
MMSEQS=mmseqs
CLOUD=./build/fb-pruner
APPEND_IDS=scripts/m8_add_ids.py

# main
BENCH_DIR=$(pwd)

# temporary directory
INSTANCE_ID=$(uuidgen)
TMP_ROOT=${BENCH_DIR}/tmp-mmseqs-plus/
TMP=${TMP_ROOT}/${INSTANCE_ID}/

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

# ==== MMSEQS CREATE DB ==== #
echo "# (1/7) Create MMSEQS target database..."
time $MMSEQS createdb 						\
	$TARGET $TARGET_MMSEQS 					\
	--shuffle 			0					\

echo "# (1/7) Create MMSEQS query database..."
time $MMSEQS createdb						\
	$QUERY $QUERY_MMSEQS 					\
	--shuffle 			0 					\

# ==== MMSEQS INDEX ======== #
echo "# (2/7) Create MMSEQS index of target database..."
KMER="7" 	# kmer length 
SPLIT="1" 	# split db in number of sub dbs
time $MMSEQS createindex 					\
	$TARGET_MMSEQS $TMP_MMSEQS				\
	-k $KMER 								\
	--split $SPLIT 							\

# ==== CLOUD INDEX ==== #
echo "# (3/7) Create CLOUD index of database..."
time $CLOUD index $TARGET $QUERY
time mv $TARGET.idx $TMP_CLOUD/target.idx 
time mv $QUERY.idx $TMP_CLOUD/query.idx

# capture database size
NUM_TARGETS=$( grep "#" -v $TMP_CLOUD/target.idx | wc -l )
NUM_QUERIES=$( grep "#" -v $TMP_CLOUD/query.idx | wc -l )
DB_SIZE=$((NUM_TARGETS))

# ==== MMSEQS SEARCH ======= $
echo "# (4/7) Performing MMSEQS Search..."
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
E_VALUE=$(( P_VALUE * DB_SIZE ))

MIN_UNGAPPED_SCORE="0"
time $MMSEQS search 
	$QUERY_MMSEQS $TARGET_MMSEQS 					\
	$RESULT_RAW_MMSEQS $TMP_MMSEQS 					\
	-k 						$KMER 					\
	--split 				$SPLIT 					\
	--min-ungapped-score 	$MIN_UNGAPPED_SCORE 	\
	-e 						$E_VALUE 				\
	--k-score 				$K_SCORE 				\
	--format-out 			$FORMAT_OUTPUT 			\
	--remove-tmp-files  	0						\

# ==== MMSEQS CONVERT == #
echo "# (5/7) Convert MMSEQS Alignments..."
time $MMSEQS convertalis							\
	$QUERY_MMSEQS $TARGET_MMSEQS  					\
	$RESULT_RAW_MMSEQS $RESULT_MMSEQS 				\

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

time $FB_PRUNER mmseqs 								\
	$TARGET_PATH $QUERY_PATH 						\
	--mmseqs-results 	$RESULT_MMSEQS				\
	--index 			$TARGET_INDEX $QUERY_INDEX 	\
	--output 			$STD_OUTPUT 				\
	--m8out 			$M8_OUTPUT 					\
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