#!/usr/bin/bash
###########################################################################
#	- FILE: 		easysearch-example.sh	
#	- DESC:  		Run easy-search pipeline examples for verify.
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

# easy search
ECHO_AND_RUN                              \
${MMORESEQS} easy-search       				  	\
	${QUERY} ${TARGET}  	   	      		  	\
  ${PREP_DIR}                   			  	\
  --mmseqs-kmer 	  6              		  	\
  --run-vitaln 	    0 								    \
	--run-postaln 	  0 								    \
	--m8out	 		      ${OUTPUT}.m8out 			\
  --myout     	    ${OUTPUT}.myout 			\
  --mythreshout	    ${OUTPUT}.mythreshout	\
	--mytimeout 	    ${OUTPUT}.mytimeout 	\
  --verbose 		    3 								    \

