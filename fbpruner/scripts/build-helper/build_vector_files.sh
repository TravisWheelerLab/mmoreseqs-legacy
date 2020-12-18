#!/bin/bash
#
# CREATES VECTOR CLASSES FROM TEMPLATE
#

# folder locations
BENCH_DIR=$(pwd)
VECTOR_DIR="../../src/objects/vectors"

# template names
UPPER_TMP="XXX"
LOWER_TMP="template"

# classes to build
ALL_UPPER=(STR INT CHAR FLT DBL RANGE TRACE BOUND)
ALL_LOWER=(str int char float double range trace bound)

cd $VECTOR_DIR

for (( J = 0; J < ${#ALLUPPER[@]}; J++ ))
do 
	LOWER=${ALL_LOWER[$J]}; 
	UPPER=${ALL_UPPER[$J]}; 
	echo "$UPPER, $LOWER"; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_$LOWER_TMP.c > vector_$LOWER.c; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_$LOWER_TMP.h > vector_$LOWER.h; 
done
