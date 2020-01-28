#!/usr/bin/bash

# Executes the Cloud Forward Backward program with test data

BENCH_DIR=$(pwd)
BUILD_DIR=$BENCH_DIR/../build/
IN_DIR=$BENCH_DIR/../data/
OUT_DIR=$BENCH_DIR/../output/

EXEC=$BUILD_DIR/cloud_fwdbck.exe
TEST_TARGET=$IN_DIR/test1_2.hmm
TEST_QUERY=$IN_DIR/test1_1.fa

$EXEC $TEST_TARGET $TEST_QUERY