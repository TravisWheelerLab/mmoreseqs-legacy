#!/usr/bin/bash
###########################################################################
#	- FILE:  prepsearch-example-1.sh	
#	- DESC:  Run step #1 of prep-search pipeline example.
###########################################################################

# input locations
BENCH_DIR=$(pwd)
MMORESEQS=${BENCH_DIR}/../build/bin/mmoreseqs
OUTPUT=${BENCH_DIR}/output/mmoreseqs.example

DB_DIR=${BENCH_DIR}/db/
QUERY=${DB_DIR}/QUERY_MSA
TARGET=${DB_DIR}/TARGET_FASTA
PREP_DIR=${BENCH_DIR}/tmp-mmoreseqs/

function ECHO_AND_RUN
{
  echo "# COMMAND: ${@}"
  echo "#"
  "${@}"
}

# prepare data files
#ECHO_AND_RUN \
${MMORESEQS} prep \
  ${QUERY} ${TARGET} \
  ${PREP_DIR} \
  --verbose 3 \
  --mmseqs-kmer 6 \
  --prep-copy 0 \

