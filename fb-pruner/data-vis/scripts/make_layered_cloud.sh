#!/bin/bash

# parse commandline
NUM_ARGS=$#
if (( $NUM_ARGS != 2 )); then
	echo "NUM_ARGS = $NUM_ARGS"
	echo "Usage: <target_dir> <query_dir>"
	exit
fi

# arguments
IN_TARGET_DIR=$1
IN_QUERY_DIR=$2

# make tmp dir 
TMP_DIR=test_output/
mkdir $TMP_DIR

# output dir
OUTPUT_DIR=final/
mkdir $OUTPUT_DIR

# program
FB_PRUNER=/home/devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/fb-pruner/build/fb-pruner

# parameters
ALPHAS=(8 10 12 14 16)
ALPHA_MAXES=(0 8 16)
BETAS=(5)

# counts
NUM_TARGET=0
NUM_QUERY=0

echo "# Beginning vizualation..."

# loop over all target/query pairs
for TARGET in $(ls $IN_TARGET_DIR)
do
	MY_TARGET=$IN_TARGET_DIR/$TARGET
	
	for QUERY in $(ls $IN_QUERY_DIR)
	do
		MY_QUERY=$IN_QUERY_DIR/$QUERY

		CUR_DIR=$OUTPUT_DIR/$TARGET/$QUERY/
		mkdir -p $CUR_DIR

		echo "target: $TARGET, query: $QUERY"
		# loop over all 
		for ALPHA in ${ALPHAS[@]}
		do 
			echo "# ALPHA = $ALPHA"
			for ALPHA_MAX in ${ALPHA_MAXES[@]}
			do
				echo "# ALPHA_MAX = $ALPHA_MAX"
				MY_ALPHA_MAX=$(($ALPHA + $ALPHA_MAX))
				for BETA in ${BETAS[@]}
				do 
					echo "# fb-pruner viz $TARGET $QUERY --alpha $ALPHA --alpha-max $MY_ALPHA_MAX --beta $BETA"
					$FB_PRUNER viz $MY_TARGET $MY_QUERY --alpha $ALPHA --alpha-max $MY_ALPHA_MAX --beta $BETA > fb-pruner-stdout.$ALPHA.$ALPHA_MAX.$BETA.txt

					# move edgebounds 
					mv $TMP_DIR/my.cloud.quad.rows.edg $CUR_DIR/cloud.$ALPHA.$ALPHA_MAX.$BETA.edg
					# move cloud mx
					mv $TMP_DIR/my.cloud_fwd.quad.mx $CUR_DIR/cloud_fwd.$ALPHA.$ALPHA_MAX.$BETA.mx 
					mv $TMP_DIR/my.cloud_bck.quad.mx $CUR_DIR/cloud_bck.$ALPHA.$ALPHA_MAX.$BETA.mx 
				done
			done
		done
	done
done

echo "# Ending vizualation."