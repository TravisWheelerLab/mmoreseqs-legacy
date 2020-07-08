#!/bin/bash
# #######################################################
#	FILE:	MAKE_ALL_BUILDS.sh
#	DESC:	Builds both DEBUG and PRODUCTION binaries
# #######################################################

mkdir bin
make clean

make BUILD=DEBUG
mv build/fb-pruner bin/fb-pruner-DEBUG 

make clean

make 
mv build/fb-pruner bin/fb-pruner-PRODUCTION