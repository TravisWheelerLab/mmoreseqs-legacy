#!/bin/bash

BENCH_DIR=$(pwd)
NUM_ARGS=$#
echo "NUM_ARGS: $NUM_ARGS"

if (( NUM_ARGS < 1 )); then
	echo "Usage: <fasta_file> <OPT=temp_file>"
	exit
fi

INPUT_FILE=$1
INPUT_FILENAME=$(basename $INPUT_FILE)
INPUT_FILEPATH=$()
OUTPUT_FILE=${INPUT_FILE}.hmm
TMP_DIR=$BENCH_DIR/tmp_convert/

if (( NUM_ARGS >= 2 )); then
	TMP_DIR=$2
fi

mkdir $TMP_DIR
rm $OUTPUT_FILE
touch $OUTPUT_FILE

NUM_MODELS=$(grep ">" $INPUT_FILE | wc -l)

GRAB_IDX=/data1/um/drich/Wheeler-Labs/benchmarks/profmark-benchmark/db/scripts/fasta_getby.py
echo  "FASTA_LENGTH: $NUM_MODELS"

for (( IDX=1; $IDX <= $NUM_MODELS; IDX=$IDX+1 ))
do
	echo "building hmm ($IDX of $NUM_MODELS)..."
	TMP_FASTA=${TMP_DIR}/${INPUT_FILENAME}.${IDX}.fa
	TMP_HMM=${TMP_DIR}/${INPUT_FILENAME}.${IDX}.hmm
	python3 $GRAB_IDX $INPUT_FILE --index $IDX --limit 1 > $TMP_FASTA
	# capture name for output file
	NAME=$(grep ">" $TMP_FASTA | cut -d' ' -f1 | cut -c2-)
	NAME_FILE=${TMP_DIR}/${NAME}
	echo "NAME: $NAME"
	mv $TMP_FASTA $NAME_FILE
	hmmbuild $TMP_HMM $NAME_FILE
	cat $TMP_HMM >> $OUTPUT_FILE
	rm $TMP_FASTA
	rm $TMP_HMM
done

rmdir $TMP_DIR
