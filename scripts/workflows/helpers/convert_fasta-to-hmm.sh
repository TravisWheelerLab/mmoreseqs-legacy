#!/bin/bash
###########################################################################
#	NAME: 		convert_fasta-to-hmm.sh	
#	AUTHOR:		David Rich
#	DESC: 		Converts file of FASTA sequences to file of single HMM profiles.
###########################################################################

# only echoes if has higher verbosity level
function echo_v
{
	if (( $VERBOSE >= $1 ))
	then
		echo $2
	fi
}
VERBOSE="${VERBOSE:-3}"

# finds the bench directory
function GET_BENCH_DIR
{
	local DIR="$(pwd)"
	echo ${DIR}
}
# get the script directory
BENCH_DIR=$(GET_BENCH_DIR)

# finds the script directory
function GET_SCRIPT_DIR
{
	local DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
	echo ${DIR}
}
SCRIPT_DIR=$(GET_SCRIPT_DIR)

# parse command line
NUM_ARGS=$#
if (( NUM_ARGS < 2 )); then
	echo "convert_fasta-to-hmm: Converts file of FASTA sequences to file of single HMM profiles."
	echo "Usage: <i:fasta_file> <o:hmm_file> | <o:tmp_folder> <hmmbuild_cmd> <splitfile_cmd>"
	exit
fi

# main file args
INPUT_FILE="$1"
OUTPUT_FILE="$2"
TMP_DIR="${3:-"${BENCH_DIR}/tmp_convert/"}"

# external programs/scripts
HMMBUILD="${4:-"hmmbuild"}"
SPLIT_FILE="${5:-"${SCRIPT_DIR}/split_fasta.py"}"

# size of input
NUM_MODELS=$(grep ">" $INPUT_FILE | wc -l)

echo_v 3 "#  INPUT_FILE:  $INPUT_FILE"
echo_v 3 "# OUTPUT_FILE:  $OUTPUT_FILE"
echo_v 3 "#  SCRIPT_DIR:  $SCRIPT_DIR"
echo_v 3 "#     TMP_DIR:  $TMP_DIR"
echo_v 3 "#  NUM_MODELS:  $NUM_MODELS"
echo_v 3 "#    HMMBUILD:  $HMMBUILD"
echo_v 3 "#  SPLIT_FILE:  $SPLIT_FILE"

mkdir ${TMP_DIR}
cp ${INPUT_FILE} ${TMP_DIR}/${INPUT_FILE}
cd ${TMP_DIR}

python ${SPLIT_FILE} ${INPUT_FILE} 1
touch ${OUTPUT_FILE}
echo_v 3 "# File split. Begin converting to hmm..."

for (( CNT = 0; CNT < NUM_MODELS; CNT += 1 ))
do
	FILE=${INPUT_FILE}.${CNT}
	echo_v 3 "#   CHILD_FILE:  $FILE"
	NAME=$(grep ">" $FILE)
	NAME=$(echo $NAME | awk '{ print $1 }' | cut -c 2- )
	echo_v 3 "#         NAME:  $NAME"
	mv ${FILE} ${NAME}
	HMM_FILE="${FILE}.hmm"
	$HMMBUILD --singlemx -o hmmbuild.tmp.out $HMM_FILE $NAME 
	cat $HMM_FILE >> $OUTPUT_FILE
	rm $NAME
	rm $HMM_FILE
done

echo_v 3 "# Files converted."
rm hmmbuild.tmp.out
rm $INPUT_FILE
mv $OUTPUT_FILE ../
cd ../
rmdir $TMP_DIR

