#!/usr/bin/sh
##################################################################
#	NAME: 	translate_ids.sh
# 	DESC:  	translate id from profile database to query database.
##################################################################

echo "in translate_ids.sh..."

#verify proper number of variables
NUM_ARGS="$#"
REQ_ARGS="5"
if (( $NUM_ARGS < $REQ_ARGS )); then 
	echo "Incorrect number of arguments."
	echo "./ <profile_db> <sequence_db> <results_db> <dummy_db> <remove_temp>"
	echo "<rm_tmp> levels: [0] keep all files [1] remove dummy database files [2] remove profile id files"
	exit
fi

PROF_DB=$1
SEQUENCE_DB=$2
RESULTS_DB=$3
DUMMY_DB=$4
RM_TEMP=$5

# script and execution location (for accessing scripts)
BENCH_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "BENCH_DIR: $BENCH_DIR"
echo "SCRIPT_DIR: $SCRIPT_DIR"

# sizes of databases
P_SIZE=$(wc -l ${PROF_DB}_h | awk '{ print $1}')
S_SIZE=$(wc -l ${SEQUENCE_DB}_h | awk '{ print $1}')
if (( P_SIZE != S_SIZE )); then
	echo "ERROR: Profile database and Sequence database are not the same size (P=$P_SIZE, Q=$S_SIZE)."
	exit
fi

echo "DUMMY_DB: $DUMMY_DB"

INDEX=${RESULTS_DB}.index
P_INDEX=${RESULTS_DB}.pid.index
S_INDEX=${RESULTS_DB}.sid.index

# create dummy alignment database ( matches profile_id to sequence_id )
python ${SCRIPT_DIR}/translate_ids/make_dummy_align_mmdb.py $DUMMY_DB $P_SIZE
# create results ( output dummy alignment database with names - this is the map )
mmseqs createtsv $PROF_DB $SEQUENCE_DB $DUMMY_DB ${DUMMY_DB}.p2s.tsv
# store profile_id file as backup
cp ${RESULTS_DB}.index ${RESULTS_DB}.pid.index
# perform translation of profile_ids to sequence_ids
python ${SCRIPT_DIR}/translate_ids/translate_mmdb_ids.py ${DUMMY_DB}.p2s.tsv ${RESULTS_DB}.pid.index ${RESULTS_DB}.sid.index
# overwrite old profile_id index file with sequence_id file
cp ${RESULTS_DB}.sid.index ${RESULTS_DB}.index

# TODO: cleanup temporary data
if (( $RM_TEMP > 0 ))
then 
	rm ${DUMMY_DB}.0
	rm ${DUMMY_DB}.index
fi

if (( $RM_TEMP > 1 ))
then 
	rm ${RESULTS_DB}.pid.index
	rm ${RESULTS_DB}.sid.index
	rm ${DUMMY_DB}.p2s.tsv
fi 
