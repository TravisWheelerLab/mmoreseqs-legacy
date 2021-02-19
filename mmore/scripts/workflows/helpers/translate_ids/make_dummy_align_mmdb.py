#!/bin/usr/python
#####################################################################
#  FILE:  make_dummydb_align.py
#  DESC:  Creates simple dummy MMseqs align database. 
# 			 Used in converting result ids from profile results to query results.
#####################################################################

# import os
import sys
# import numpy as np

#####################################################################
##  MAIN  ###########################################################
#####################################################################

if len(sys.argv) < 3:
	print("./ <db_name> <i:db_size>")
	exit(1)

my_args = {}
my_args["db_name"] = sys.argv[1]
my_args["db_size"] = int(sys.argv[2])

my_args["db_index"] = my_args["db_name"] + ".index"
my_args["db_0"] = my_args["db_name"] + ".0"

fp_index = open(my_args["db_index"], "w")
fp_0 = open(my_args["db_0"], "w")

total_offset = 0
entry_size = 0

for i in range(my_args["db_size"]):
	line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\x00".format(i, i, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	entry_size = len(line)
	fp_0.write(line)

	index = "{} {} {}\n".format(i, total_offset, entry_size)
	total_offset += entry_size
	fp_index.write(index)

fp_0.close()
fp_index.close()

print("# completed successfully.")