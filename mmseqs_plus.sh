#!/usr/local/bin/bash
#
#	NAME: 		mmseqs_plus.sh	
#	AUTHOR:		David Rich
#	DESC: 		runs mmseqs, then pipes output to fb-pruner search
#

BENCHDIR=$(pwd)

mkdir $(BENCHDIR)/tmp-mmseqs-plus/

mmseqs 