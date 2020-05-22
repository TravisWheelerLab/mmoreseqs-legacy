#!/usr/bin/bash

# Executes the HMMER program with test data

BENCH_DIR=$(pwd)

HMMER_DIR=$BENCH_DIR/../../tools/hmmer_test/
cd $HMMER_DIR
echo $(pwd)

BUILD_DIR=/bin/
IN_DIR=/tutorial/
OUT_DIR=/output/

EXEC=./$BUILD_DIR/hmmsearch

# Segfault error?
# TEST_TARGET=$IN_DIR/test1_2.hmm
# TEST_QUERY=$IN_DIR/test1_1.fa
TEST_TARGET=$IN_DIR/test_2.hmm 
TEST_QUERY=$IN_DIR/test_1.fa

head $TEST_TARGET
head $TEST_QUERY

echo $EXEC --cpu 0 $TEST_TARGET $TEST_QUERY
$EXEC --cpu 0 $TEST_TARGET $TEST_QUERY

# cd $BENCH_DIR