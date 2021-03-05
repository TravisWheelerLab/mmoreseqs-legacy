#!/bin/usr/env python
#####################################################################
#  FILE:  convert_result_ids
#  DESC:  Convert result ids from profile results to query results.
#####################################################################

import os,sys
import numpy as np

def convert_prefilter_results():
	return

def convert_align_results():
	return

#####################################################################
##  MAIN  ###########################################################
#####################################################################

my_args = {}

my_ids = {}
my_ids["profile"] = {}
my_ids["query"] = {}
my_ids["p2q"] = {}
my_ids["q2p"] = {}
for key in my_ids.keys():
	my_ids[key]["n2i"] = {}
	my_ids[key]["i2n"] = {}

req_args = 4
if len(sys.argv) < req_args:
	print("# ERROR: incorrect number of args.")
	print("# ./ <i:lookup_table> <i:db_index> <o:db_index>")
	exit(-1)

# for i in range(len(sys.argv)):
# 	print("# {}: {}".format(i, sys.argv[i]))

# load args
my_args["lookup_table"] = sys.argv[1]
my_args["db_index_in"] = sys.argv[2]
my_args["db_index_out"] = sys.argv[3]

# load lookup table
fp_lookup = open(my_args["lookup_table"], "r")
line = fp_lookup.readline()
while line:
	fields = line.split()
	p_name = fields[0]
	q_name = fields[1]
	my_id  = fields[2]

	my_ids["profile"]["n2i"][p_name] = my_id
	my_ids["query"]["n2i"][q_name] = my_id
	my_ids["profile"]["n2i"][my_id] = p_name
	my_ids["query"]["n2i"][my_id] = q_name

	line = fp_lookup.readline()
fp_lookup.close()

# build id translation table
names = my_ids["profile"]["n2i"].keys()
for i in range(len(names)):
	name = names[i]
	p_id = my_ids["profile"]["n2i"][name]
	q_id = my_ids["query"]["n2i"][name]

	my_ids["p2q"][p_id] = q_id
	my_ids["q2p"][q_id] = p_id


# translate ids in database index
fp_index_in = open(my_args["db_index_in"], "r")
fp_index_out = open(my_args["db_index_out"], "w")

in_line = fp_index_in.readline()
while in_line:
	#in_line = in_line.strip()
	#in_line = in_line.replace('\x00', '')
	fields = in_line.split()
	if len(fields) < 3:
		print(r"# ERROR: only {} fields. LINE = {}".format(len(fields), in_line))
		in_line = fp_index_in.readline()
		continue
		
	p_id = fields[0]
	offset = fields[1]
	size = fields[2]

	q_id = my_ids["p2q"][p_id]
	out_line = "{}\t{}\t{}\n".format(q_id, offset, size)
	fp_index_out.write(out_line)
	
	in_line = fp_index_in.readline()

print("# completed successfully.")
