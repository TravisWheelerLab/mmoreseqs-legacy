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
# (1) mmseqs creates index of query/target files
# (2) fb-pruner creates index of query/target files
# (3) mmseqs search
# (4) 

# programs
MMSEQS=mmseqs
CLOUD=fb-pruner

# main directories
BENCH_DIR=$(pwd)
INSTANCE_ID=$(uuidgen)
TMP_ROOT=${BENCH_DIR}/tmp-mmseqs-plus/
TMP=${TMP_ROOT}/${INSTANCE_ID}/

# result destination files
MMSEQS_RES=${TMP}/mmseqs.m8
CLOUD_RES=${TMP}/cloud.m8

# bash variables (verify proper number of variables)
if [$# -ne 1]; then 
	echo "illegal number of parameters"
	exit(1)
fi
QUERY_PATH=$1
TARGET_PATH=$2
RESULTS=$3

# temporary example
QUERY_PATH=${BENCH_DIR}/test_input/3-PAP.fa
TARGET_PATH=${BENCH_DIR}/test_input/3-PAP.hmm

# make tmp directories
mkdir $TMP_ROOT
mkdir $TMP
mkdir $CLOUD_TMP
mkdir $MMSEQS_TMP

# ==== MMSEQS INDEX ==== #

# ==== MMSEQS SEARCH ==== #

# mmseqs search
# NOTE: need to reveal ids
# --format-output STR       	CSV list of output columns from: query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid,taxid,taxname,taxlineage  [query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits]
# -s FLOAT                     	Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [5.700]
# -k INT                       	k-mer length (0: automatically set to optimum) [0]
# --k-score INT                	k-mer threshold for generating similar k-mer lists [2147483647]
# --max-rejected INT           	Maximum rejected alignments before alignment calculation for a query is stopped [2147483647]
# --max-accept INT             	Maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
#  --local-tmp STR              Path where some of the temporary files will be created []
# --remove-tmp-files BOOL      Delete temporary files [1]
# --e-profile FLOAT            Include sequences matches with < e-value thr. into the profile (>=0.0) [0.001]

FORMAT_OUTPUT="--format-output qsetid,tsetid,qset,tset,query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
SENSITIVITY="--s 7.5"
TMP_DEL="--remove-tmp-files 0"
EVAL_THRESHOLD="--e-profile 10000.0"

echo "(1/2) PERFORMING MMSEQS SEARCH..."
time $MMSEQS easy-search $QUERY_PATH $DB_PATH $MMSEQS_RESULTS $MMSEQS_TMP $MMSEQS $FORMAT_OUTPUT $SENSITIVITY $TMP_DEL $EVAL_THRESHOLD
echo "(1/2) ...MMSEQS SEARCH COMPLETED."

# ==== CLOUD INDEX ==== #

CLOUD_TARGET_INDEX=$CLOUD_TMP/target.idx
CLOUD_QUERY_INDEX=$CLOUD_TMP/query.idx
CLOUD_INDEX_IN="--index $CLOUD_TARGET_INDEX $CLOUD_QUERY_INDEX"

echo "(3) PERFORMING CLOUD INDEXING..."
time $CLOUD index $TARGET_PATH $QUERY_PATH $CLOUD_INDEX_IN
echo "(3) ...CLOUD INDEXING COMPLETED."

# cloud search options:
# --alpha FLOAT 				alpha cloud tuning parameter
# --beta INT 					beta cloud tuning parameter
# --pipeline STR/INT 			Choose: [test,main,mmseqs,index]
# --threshold FLOAT 			Score Reporting Threshold
# --format-input STR 			CSV list of input format
# --mmseqs-results STR 			Location of mmseqs results
# --mmseqs-tid STR 				Location of mmseqs target id lookup table
# --mmseqs-qid STR 				Location of mmseqs query id lookup table

MMSEQS_RES_IN="--mmseqs-results $MMSEQS_RESULTS" 
MMSEQS_LOOKUP_IN="--mmseqs-lookup $MMSEQS_LOOKUP"
CLOUD_INDEX_IN="--index-input $CLOUD_INDEX"
RESULTS_OUT="--output $CLOUD_RES"

echo "(4) PERFORMING CLOUD SEARCH..."
$CLOUD mmseqs $TARGET_PATH $QUERY_PATH $MMSEQ
echo "(4) ...CLOUD SEARCH COMPLETED."

# remove temporary directories