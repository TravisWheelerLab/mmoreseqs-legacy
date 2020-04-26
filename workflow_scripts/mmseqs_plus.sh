 #!/usr/local/bin/bash
 ###############################################################################
 #     FILE:	mmseqs_plus.sh
 #  PURPOSE:	run mmseqs plus
 #
 #   AUTHOR: 	David Rich
 #      BUG:
 ###############################################################################

# programs 
MMSEQS_SEARCH=mmseqs
CLOUD_SEARCH=fb-pruner
MMSEQS_PLUS="fb-pruner mmseqs"
# directories
BENCH_DIR=$(pwd)/
TMP_ROOT=$BENCH_DIR/tmp/
TMP=$BENCH_DIR/$TMP_ROOT/latest/
# input files
TARGET=$BENCH_DIR/db/pmark.hmm
QUERY=$BENCH_DIR/db/pmark.fa
T_INDEX=$BENCH_DIR/db/pmark.hmm.idx
Q_INDEX=$BENCH_DIR/db/pmark.fa.idx
RESULTS_IN=$BENCH_DIR/results/mmseqs_res.m8+
# output files
OUTFILE=$BENCH_DIR
# parameters
ALPHA=10.0
BETA=5
# range of outfile lines to read in
TOTAL_RESULTS=465779
RANGE_BEG=0
RANGE_END=5

# commandline args
NUM_JOBS=100
BATCH_SIZE=(($TOTAL_RESULTS/))
ID=$1
RANGE_BEG=(())


# mmseqs plus call
echo "running mmseqs plus..."
echo "$MMSEQS_PLUS $TARGET $QUERY --index $TARGET $QUERY --mmseqs-m8+ $RESULTS_IN --alpha $ALPHA --beta $BETA --mmseqs-range $RANGE_BEG $RANGE_END"
$MMSEQS_PLUS $TARGET $QUERY --index $TARGET $QUERY --mmseqs-m8+ $RESULTS_IN --alpha $ALPHA --beta $BETA --mmseqs-range $RANGE_BEG $RANGE_END
