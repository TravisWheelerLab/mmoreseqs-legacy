#!/usr/bin/bash
###########################################################################
# - FILE:  prepsearch-example.sh	
# - DESC:  Run prep-search pipeline examples.
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

# search based on prepared files
#ECHO_AND_RUN \
${MMORESEQS} prep-search \
  ${PREP_DIR} \
  --verbose 3 \
  \
  --mmseqs-kmer 6 \
  --run-vitaln 0 \
  --run-postaln 1 \
  --m8out ${OUTPUT}.m8out \
  --myout ${OUTPUT}.myout \
  --mythreshout ${OUTPUT}.mythreshout \
