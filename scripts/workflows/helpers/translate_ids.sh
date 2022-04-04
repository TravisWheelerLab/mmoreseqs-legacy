#!/usr/bin/sh
##################################################################
#	 - NAME: 	translate_ids.sh
#  - DESC:  Translate id from profile database to query database.
##################################################################

#verify proper number of variables
NUM_ARGS="$#"
REQ_ARGS="5"
if (( $NUM_ARGS < $REQ_ARGS )); then 
	echo "Incorrect number of arguments ($NUM_ARGS of $REQ_ARGS)."
	echo "./translate_ids.sh <profile_db> <sequence_db> <results_db> <dummy_db> <remove_temp>"
	echo "<rm_tmp> levels: [0] keep all files [1] remove dummy database files [2] remove profile id files"
	echo "./translate_ids.sh [1]$1 [2]$2 [3]$3 [4]$4 [5]$5"
	exit
fi

PROF_DB=$1
SEQ_DB=$2
RESULTS_DB=$3
DUMMY_DB=$4
IDENTITY_DB=$5
RM_TEMP="${6:-0}"
VERBOSE="${7:-3}"

if (( $VERBOSE >= 3 ))
then 
{
	echo "PROF_DB: $PROF_DB"
	echo "SEQ_DB: $SEQ_DB"
	echo "RESULTS_DB: $RESULTS_DB"
	echo "DUMMY_DB: $DUMMY_DB"
	echo "IDENTITY_DB: $IDENTITY_DB"
	echo "RM_TEMP: $RM_TEMP"
}
fi

# script and execution location (for accessing scripts)
BENCH_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# sizes of databases
P_SIZE=$(wc -l ${PROF_DB}_h | awk '{ print $1}')
S_SIZE=$(wc -l ${SEQ_DB}_h | awk '{ print $1}')
if (( P_SIZE != S_SIZE )); then
	echo "ERROR: Profile database and Sequence database are not the same size (P=$P_SIZE, Q=$S_SIZE)."
	exit
fi

INDEX="${RESULTS_DB}.index"
P_INDEX="${RESULTS_DB}.pid.index"
S_INDEX="${RESULTS_DB}.sid.index"

# create dummy alignment database (matches profile_id to sequence_id by id value)
python ${SCRIPT_DIR}/translate_ids/make_dummy_align_mmdb.py $DUMMY_DB $P_SIZE
# generate results ( output dummy alignment database with names... this is the map of id->name, id is line number )
mmseqs createtsv $PROF_DB $SEQ_DB $DUMMY_DB ${DUMMY_DB}.n2i.mm_tsv --threads 1
mmseqs convertalis $PROF_DB $SEQ_DB $DUMMY_DB ${DUMMY_DB}.mm_m8 --threads 1 -v 3
# store profile_id file as backup
cp ${RESULTS_DB}.index ${RESULTS_DB}.pid.index
# perform translation of profile_ids to sequence_ids
python ${SCRIPT_DIR}/translate_ids/translate_mmdb_ids.py 	\
${DUMMY_DB}.mm_m8 			${DUMMY_DB}.i2i.mm_tsv 				\
${RESULTS_DB}.pid.index 	${RESULTS_DB}.sid.index				\
$IDENTITY_DB
# generate results ( output dummy alignment database with names... this is the map of pid->sid )
# mmseqs createtsv $PROF_DB $SEQ_DB $DUMMY_DB ${DUMMY_DB}.n2n.mm_tsv --threads 1
mmseqs convertalis $PROF_DB $SEQ_DB $IDENTITY_DB ${IDENTITY_DB}.mm_m8 --threads 1 -v 3
# # overwrite old profile_id index file with sequence_id file
# cp ${RESULTS_DB}.sid.index ${RESULTS_DB}.index
# create alignment database which matches profiles to sequences
# python ${SCRIPT_DIR}/translate_ids/make_identity_align_mmdb.py $IDENTITY_DB ${DUMMY_DB}.i2i.mm_tsv

# TODO: cleanup temporary data
if (( $RM_TEMP > 0 ))
then 
	rm ${DUMMY_DB}
	rm ${DUMMY_DB}.0
	rm ${DUMMY_DB}.index
	rm ${DUMMY_DB}.dbtype
fi

if (( $RM_TEMP > 1 ))
then 
	# rm ${RESULTS_DB}.pid.index
	# rm ${RESULTS_DB}.sid.index
	rm ${DUMMY_DB}.p2s.tsv
fi 
