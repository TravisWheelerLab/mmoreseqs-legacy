#!/bin/bash

BENCH_DIR=$(pwd)
NUM_ARGS=$#
echo "   NUM_ARGS:  $NUM_ARGS"

if (( NUM_ARGS < 1 )); then
	echo "Usage: <fasta_file> <OPT=temp_file>"
	exit
fi

# get location of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# external programs/scripts
SPLIT_FILE=$SCRIPT_DIR/split_fasta.py
HMMBUILD=hmmbuild

INPUT_FILE=$1
INPUT_FILENAME=$(basename $INPUT_FILE)
INPUT_FILEPATH=$()
OUTPUT_FILE=${INPUT_FILE}.hmm
TMP_DIR=$BENCH_DIR/tmp_convert/

NUM_MODELS=$(grep ">" $INPUT_FILE | wc -l)

echo "#  INPUT_FILE:  $INPUT_FILE"
echo "# OUTPUT_FILE:  $OUTPUT_FILE"
echo "#  SCRIPT_DIR:  $SCRIPT_DIR"
echo "#     TMP_DIR:  $TMP_DIR"
echo "#  NUM_MODELS:  $NUM_MODELS"

mkdir $TMP_DIR
cp $INPUT_FILE $TMP_DIR/$INPUT_FILE
cd $TMP_DIR

python $SPLIT_FILE $INPUT_FILE 1
touch $OUTPUT_FILE
echo "# File split. Begin converting to hmm..."

for (( CNT = 0; CNT < NUM_MODELS; CNT += 1 ))
do
	FILE=${INPUT_FILE}.${CNT}
	echo "#   CHILD_FILE:  $FILE"
	NAME=$(grep ">" $FILE)
	NAME=$(echo $NAME | awk '{ print $1 }' | cut -c 2- )
	echo "#         NAME:  $NAME"
	mv $FILE $NAME
	HMM_FILE=$FILE.hmm
	hmmbuild --singlemx -o hmmbuild.tmp.out $HMM_FILE $NAME 
	cat $HMM_FILE >> $OUTPUT_FILE
	rm $NAME
	rm $HMM_FILE
done

echo "# Files converted."
rm hmmbuild.tmp.out
rm $INPUT_FILE
mv $OUTPUT_FILE ../
cd ../
rmdir $TMP_DIR

