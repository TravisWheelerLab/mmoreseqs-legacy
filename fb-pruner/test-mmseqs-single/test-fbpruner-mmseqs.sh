#!/usr/bin/bash
###############################################################################
#    FILE:  mmseqs-plus-test
#   BRIEF:  Test fbpruner's mmseqs-plus pipeline
#  AUTHOR:  Dave Rich
###############################################################################

QUERY=test1_2.hmm
TARGET=test1_1.fa

QUERY_IDX=test1_2.hmm.idx
TARGET_IDX=test1_1.fa.idx

MMSEQS_M8=mmseqs_n1.test.n1.m8

ALPHA=6 
BETA=8
GAMMA=5

OUTPUT=test-results.out
TBLOUT=test-results.tblout
M8OUT=test-results.m8
MYOUTPUT=test-results.myout

VERBOSE=1

time ../build/fb-pruner mmseqs 					\
					$QUERY $TARGET 				\
	--index 		$QUERY_IDX $TARGET_IDX 		\
	--alpha 		$ALPHA 						\
	--beta 			$BETA						\
	--gamma 		$GAMMA 						\
	--mmseqs-m8 	$MMSEQS_M8 					\
	--m8out 		$M8OUT 						\
 	--verbose 		$VERBOSE					\
												\

