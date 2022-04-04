#!/bin/usr/python
#####################################################################
#  - FILE:  make_line_index.py
#  - DESC:   Creates a simple line index of file.
# 			 Series of line entries {line_length} {total_offset}
#####################################################################

# import os
import sys
# import numpy as np

# identifiers for different database types
my_dbtypes = {}
my_dbtypes["amino"] 			= ( "%c%c%c%c" % (0,0,0,0) )
my_dbtypes["nucleotide"] 	= ( "%c%c%c%c" % (1,0,0,0) )
my_dbtypes["profile"] 		= ( "%c%c%c%c" % (2,0,0,0) )
my_dbtypes["prefilter"] 	= ( "%c%c%c%c" % (7,0,0,0) )
my_dbtypes["align"] 			= ( "%c%c%c%c" % (5,0,0,0) )
my_dbtypes["generic"] 		= ( "%c%c%c%c" % (12,0,0,0) )

#####################################################################
##  MAIN  ###########################################################
#####################################################################

if len(sys.argv) < 3:
	print("./ <i:db_file> <o:index_file>")
	exit(1)

my_args = {}
my_args["db_type"] = "align"

my_args["db_file"] 		= sys.argv[1]
my_args["index_file"] 	= sys.argv[2]

my_files = {}
my_files["db_file"] 		= "{}".format(my_args["db_file"])
my_files["index_file"] 	= "{}".format(my_args["index_file"])

fp_db 	= open(my_files["db_file"], "r")
fp_index	= open(my_files["index_file"], "w")

total_offset = 0
entry_size = 0

for line in fp_db:
	# get next line
	sline = line.strip()
	entry_size = len(line)

	# make entry for index
	out_line = "{}\t{}\n".format(total_offset+1, entry_size)
	total_offset += entry_size
	fp_index.write(out_line)

fp_db.close()
fp_index.close()

print("# completed successfully.")
