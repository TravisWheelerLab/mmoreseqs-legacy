#!/bin/bash

# working directory
BENCH_DIR=$(pwd) 

# location of data-vis program
DATA_VIS=${BENCH_DIR}/../data-vis/scripts/build_matrix_layered_edgebounds.py
DEFAULT_OUTPUT=layer-map.LAYMAP.jpg

# query and target lists
QUERIES=(0 1 2 3 4 5)
TARGETS=(0 1 2 3 4 5)
ALPHAS=(8 10 12 14 16)
ALPHA_MAXES=(0 8 16)

for Q in ${QUERIES[@]}
do
	for T in ${TARGETS[@]}
	do
		cd final/target.1k.${T}.hmm/query.1k.${Q}.fasta/

		for A in ${ALPHAS[@]}
		do
			# create alpha slices
			python $DATA_VIS --no-show $( ls cloud.${A}.*.5.edg )
			mv $DEFAULT_OUTPUT $BENCH_DIR/layers.${Q}.${T}.alpha.${A}.jpg
		done

		for M in ${ALPHA_MAXES[@]}
		do	
			# create alpha-max slice
			python $DATA_VIS --no-show $(ls cloud.*.${M}.5.edg )
			mv $DEFAULT_OUTPUT $BENCH_DIR/layers.${Q}.${T}.alphamax.${M}.jpg
		done

		cd $BENCH_DIR
	done
done


