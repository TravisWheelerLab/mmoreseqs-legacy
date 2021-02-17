#!/usr/bin/sh
######################################################
#	NAME: 	translate_ids.sh
# 	DESC:  	translate id from profile database to query database.
######################################################

REQ_ARGS=3
NUM_ARGS=$#

if (( NUM_ARGS < REQ_ARGS )); then
	echo "Incorrect number of arguments."
	echo "./ <profile_db> <query_db> <results_db>"
	exit
fi

PROF_DB=$1
QUERY_DB=$2
RESULTS_DB=$3
TEST_DB="testdb"

BENCH_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "SCRIPT DIR: $SCRIPT_DIR"

P_SIZE=$(wc -l ${PROF_DB}_h | awk '{ print $1}')
Q_SIZE=$(wc -l ${QUERY_DB}_h| awk '{ print $1}')
echo "P_SIZE: $P_SIZE, Q_SIZE: $Q_SIZE"
if (( P_SIZE != Q_SIZE )); then
	echo "ERROR: Profile database and Query database are not the same size (P=$P_SIZE, Q=$Q_SIZE)."
	exit
fi

python $SCRIPT_DIR/make_testdb_align.py $TEST_DB $P_SIZE

mmseqs createtsv $PROF_DB $QUERY_DB $TEST_DB ${TEST_DB}.p2q.tsv

if [ ! -f ${RESULTS_DB}.profile.index ]; then
	cp ${RESULTS_DB}.index ${RESULTS_DB}.profile.index
fi

python $SCRIPT_DIR/translate_result_ids.py ${TEST_DB}.p2q.tsv ${RESULTS_DB}.profile.index ${RESULTS_DB}.query.index

cp ${RESULTS_DB}.query.index ${RESULTS_DB}.index

# TODO: cleanup temporary data
#rm testdb.0
#rm testdb.index

