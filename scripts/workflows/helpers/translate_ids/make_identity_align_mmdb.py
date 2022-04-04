#!/bin/usr/python
#####################################################################
#  - FILE:  make_identity_align_mmdb.py
#  - DESC:   Creates simple dummy MMseqs align database. 
# 			 Used in converting result ids from profile results to query results.
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
	print("./ <o:db_name> <i:id2id_tsv> <i:dbtype>")
	exit(1)

my_args = {}
my_args["db_type"] 	= "align"

my_args["db_name"] 	= sys.argv[1]
my_args["id_tsv"] 	= sys.argv[2]
if len(sys.argv) > 3:
	my_args["db_type"] 	= sys.argv[3]

my_files = {}
my_files["id_tsv"]  	= "{}".format(my_args["id_tsv"])
my_files["db_index"] = "{}.index".format(my_args["db_name"])
my_files["db_0"] 		= "{}".format(my_args["db_name"])
my_files["db_type"] 	= "{}.dbtype".format(my_args["db_name"])

fp_id 	= open(my_files["id_tsv"], "r")

fp_index = open(my_files["db_index"], "w")
fp_0 		= open(my_files["db_0"], "w")
fp_type 	= open(my_files["db_type"], "w")

total_offset 	= 0
entry_size 		= 0

id_line = fp_id.readline()

# build a db entry for every entry in id2id tsv file
while id_line:
	id_fields 	= id_line.split()
	p_id 		 	= id_fields[0]
	q_id 		 	= id_fields[1]
	i 		 	 	= id_fields[2]
	name 			= id_fields[3]

	# add entry for database (q_id is main key)
	# NOTE: do I need a EOF (x00) at end of each line?
	line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\x00".format(q_id, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	entry_size = len(line)
	fp_0.write(line)

	# add entry for the index (p_id is main key)
	index = "{}\t{}\t{}\n".format(p_id, total_offset, entry_size)
	total_offset += entry_size
	fp_index.write(index)

	id_line = fp_id.readline()

# make dbtype file
fp_type.write(my_dbtypes[my_args["db_type"]])

fp_0.close()
fp_index.close()
fp_type.close()

print("# completed successfully.")
