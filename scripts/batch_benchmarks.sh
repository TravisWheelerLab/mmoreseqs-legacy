 #!/usr/local/bin/bash

 ###############################################################################
 #  @file  batch benchmarks
 #  @brief Render benchmark scores as a line graph
 #
 #  @author Dave Rich
 #  @bug Lots.
 ###############################################################################


BATCH_BENCHMARK() 
{
   CLOUD_SEARCH="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/build/cloud_fwdbck.exe"
   TEST_DIR="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data/test_list/"
   STATS_DIR="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/stats/"
   DATAVIS_DIR="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data-vis/alpha-linegraphs/"
   DOMTBLOUT="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/phmmer_001.domtblout"

   echo "HMM_FILE: $HMM_FILE, FA_FILE: $FA_FILE"

   # rm $STATS
   # touch $STATS

   for (( ALPH = 2; ALPH < 20; ALPH += 2 ))
   do
      ${CLOUD_SEARCH} ${TEST_DIR}/${HMM_FILE} ${TEST_DIR}/${FA_FILE} -o ${STATS_DIR}/${FA_FILE} -a ${ALPH}
   done

   BUILD_LINEGRAPH="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/scripts/build_benchmark_linegraph.py"

   python ${BUILD_LINEGRAPH} ${STATS_DIR}/${FA_FILE}
}

if (( $# != 2 )); 
then 
{
   echo "Usage: <hmm_file> <fasta_file>"
   exit 1
} 
fi

HMM_FILE=$1
FA_FILE=$2

BATCH_BENCHMARK