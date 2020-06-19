#!/bin/bash
##################################################################################
#	NAME: 	MMSEQS_REBUILD.sh
# 	DESC: 	Script for performing a complete rebuild of the mmseqs binaries
#		from source.
# 	AUTH: 	David Rich
##################################################################################
set -x #echo on


rm -r build
#git clone https://github.com/soedinglab/MMseqs2.git
#cd MMseqs2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make -j 4
make install
export PATH=$(pwd)/bin/:$PATH
