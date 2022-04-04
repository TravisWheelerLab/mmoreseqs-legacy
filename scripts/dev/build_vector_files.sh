#!/bin/bash
##############################################################
# - FILE: build_vector_files.sh
# Desc: Creates Vector files from template.
##############################################################

echo "Running 'build_vector_files.sh'..."

# folder locations
BENCH_DIR=$(pwd)
VECTOR_DIR="../../src/objects/vectors"

# template names
UPPER_TMP="XXX"
LOWER_TMP="template"

# classes to build
ALL_UPPER=(STR INT PTR CHAR FLT DBL RANGE TRACE BOUND)
ALL_LOWER=(str int ptr char float double range trace bound)

cd $VECTOR_DIR

LEN=${#ALL_UPPER[@]}
for (( J = 0; J < $LEN; J++ ))
do
	LOWER=${ALL_LOWER[$J]}; 
	UPPER=${ALL_UPPER[$J]}; 

	echo "building 'vector_${LOWER}.c' from '__vector_${LOWER_TMP}__.c'..."	
	echo "$UPPER, $LOWER"; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_${LOWER_TMP}.c > vector_${LOWER}.c; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_${LOWER_TMP}.h > vector_${LOWER}.h; 
done
