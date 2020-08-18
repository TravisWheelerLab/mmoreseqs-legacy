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
ALL_UPPER=(INT CHAR FLT DBL RANGE TRACE BOUND XXX)
ALL_LOWER=(int char float double range trace bound template)

cd $VECTOR_DIR

for (( J=0; J<7; J++ ))
do 
	LOWER=${ALL_LOWER[$J]}; 
	UPPER=${ALL_UPPER[$J]}; 
	echo "$UPPER, $LOWER"; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_$LOWER_TMP.c > vector_$LOWER.c; 
	sed "s/$UPPER_TMP/$UPPER/g; s/vector_$LOWER_TMP/vector_$LOWER/g" vector_$LOWER_TMP.h > vector_$LOWER.h; 
done
