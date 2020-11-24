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

VERBOSE=1

time ../../more-tools/hmmer-test/bin/hmmsearch 	\
	--cpu 			0 							\
	--max 										\
					$QUERY $TARGET 				\
												\

