 #!/usr/local/bin/bash

BATCH="bash /Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/scripts/batch_benchmarks.sh"
TEST_DIR="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data/test_list"
BENCH_DIR="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner"

PWD=$(pwd)

cd $TEST_DIR
HMM_FILES=$(ls *.hmm)

for HMM in $HMM_FILES
do
   FNAME=$(basename -- "$HMM")
   EXT="${FNAME##*.}"
   FNAME="${FNAME%.*}"
   FA_FILES=$(ls ${FNAME}.*.fa)

   cd $BENCH_DIR

   for FA in $FA_FILES
   do 
      echo "$BATCH $HMM $FA"

      # TODO:check if file already exists!!

      $BATCH $HMM $FA
   done

   cd $TEST_DIR

done