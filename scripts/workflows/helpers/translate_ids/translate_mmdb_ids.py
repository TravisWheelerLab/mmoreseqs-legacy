#!/bin/usr/env python
#####################################################################
#  - FILE:  translate_mmdb_ids.py
#  - DESC:   Trnslate mmdb from profile to query ids.
#####################################################################

from __future__ import print_function
import os,sys
import numpy as np

# identifiers for different database types
my_dbtypes = {}
my_dbtypes["amino"] 			= ( "%c%c%c%c" % (0,0,0,0) )
my_dbtypes["nucleotide"] 	= ( "%c%c%c%c" % (1,0,0,0) )
my_dbtypes["profile"] 		= ( "%c%c%c%c" % (2,0,0,0) )
my_dbtypes["prefilter"] 	= ( "%c%c%c%c" % (7,0,0,0) )
my_dbtypes["align"] 			= ( "%c%c%c%c" % (5,0,0,0) )
my_dbtypes["generic"] 		= ( "%c%c%c%c" % (12,0,0,0) )

# error printing 
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# load lookup table
def load_name_lookup_table(my_files, my_ids):
	fp_lookup = open(my_files["name_tbl_in"], "r")
	line = fp_lookup.readline()
	line_number = 0
	# get file line by line
	while line:
		# parse line
		fields = line.split()
		p_name = fields[0]
		q_name = fields[1]
		# my_id  = fields[2]
		my_id  = line_number

		# load name->id and id->name dictionaries
		my_ids["profile"]["n2i"][p_name] = my_id
		my_ids["query"]["n2i"][q_name] 	= my_id
		my_ids["profile"]["i2n"][my_id] 	= p_name
		my_ids["query"]["i2n"][my_id] 	= q_name

		line = fp_lookup.readline()
		line_number += 1
	fp_lookup.close()
	return

# build query->profile translation table
def build_translation_table(my_files, my_ids):
	
	fp_idlookup = open(my_files["id_tbl_out"], "w")
	names = list(my_ids["profile"]["n2i"].keys())

	for i in range(len(names)):
		
		name = names[i]
		p_id = my_ids["profile"]["n2i"][name]
		q_id = my_ids["query"]["n2i"][name]

		my_ids["p2q"][p_id] = q_id
		my_ids["q2p"][q_id] = p_id
		
		line = "{}\t{}\t{}\t{}".format(p_id, q_id, i, name)
		fp_idlookup.write( "{}\n".format(line) )
		eprint("[{}] {}: {} -> {}".format(i, name, p_id, q_id) )

	fp_idlookup.close()
	return

# translate db_index (profile -> query)
def translate_db_index(my_files, my_ids):
	# translate ids in database index from profile ID -> query/sequence ID
	fp_index_in 	= open(my_files["db_index_in"], "r")
	fp_index_out 	= open(my_files["db_index_out"], "w")

	in_line = fp_index_in.readline()
	while in_line:
		#in_line = in_line.strip()
		#in_line = in_line.replace('\x00', '')
		fields = in_line.split()
		if len(fields) < 3:
			eprint(r"# ERROR: only {} fields. LINE = {}".format(len(fields), in_line))
			in_line = fp_index_in.readline()
			continue
			
		p_id 		= int(fields[0])
		offset 	= int(fields[1])
		size 		= int(fields[2])

		q_id = my_ids["p2q"][p_id]
		out_line = "{}\t{}\t{}\n".format(q_id, offset, size)
		fp_index_out.write(out_line)
		
		in_line = fp_index_in.readline()
	return 

# translate db_data (query -> profile)
def translate_db_data(my_files, my_ids):
	# translate ids in database index from profile ID -> query/sequence ID
	fp_index_in 	= open(my_files["db_index_in"], "r")
	fp_index_out 	= open(my_files["db_index_out"], "w")

	in_line = fp_index_in.readline()
	while in_line:
		#in_line = in_line.strip()
		#in_line = in_line.replace('\x00', '')
		fields = in_line.split()
		if len(fields) < 3:
			eprint(r"# ERROR: only {} fields. LINE = {}".format(len(fields), in_line))
			in_line = fp_index_in.readline()
			continue
			
		q_id 	= int(fields[0])
		p_id 	= my_ids["q2p"][q_id]

		out_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\x00".format(
			p_id,      fields[1], fields[2], fields[3], fields[4], 
			fields[5], fields[6], fields[7], fields[8], fields[9])
		fp_index_out.write(out_line)
		
		in_line = fp_index_in.readline()
	return 

# make identity database
def make_identity_database(my_files, my_ids):
	fp_data 	= open(my_files["db_identity_data"], "w")
	fp_index = open(my_files["db_identity_index"], "w")
	fp_type 	= open(my_files["db_identity_type"], "w")

	total_offset 	= 0
	entry_size 		= 0

	p_keys 	= set(my_ids["profile"]["n2i"].keys())
	p2q_keys = list(my_ids["p2q"].keys())
	eprint(p2q_keys)

	for i in range(len(p_keys)):
		q_id = my_ids["p2q"][i]

		# make entry for database
		# targetID  alnScore  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen
		# [queryOrfStart] [queryOrfEnd] [dbOrfStart] [dbOrfEnd] [alnCigar]
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\x00".format(q_id, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
		entry_size = len(line)
		fp_data.write(line)

		# make entry for index
		index = "{}\t{}\t{}\n".format(i, total_offset, entry_size)
		total_offset += entry_size
		fp_index.write(index)

	# make dbtype file
	fp_type.write( my_dbtypes["align"] )	

	return


#####################################################################
##  MAIN  ###########################################################
#####################################################################

my_ids = {}
my_ids["profile"] = {}
my_ids["query"] 	= {}
my_ids["p2q"] 		= {}
my_ids["q2p"] 		= {}

for key in my_ids.keys():
	my_ids[key]["n2i"] = {}
	my_ids[key]["i2n"] = {}

req_args = 4
if len(sys.argv) < req_args:
	print("# ERROR: incorrect number of args.")
	print("# ./ <i:name_lookup_tbl> <i:db_index> <i:db_data> <o:id_lookup_table> <o:db_index> <o:db_data>")
	exit(1)

# load args
my_args = {}
my_args["name_tbl_in"] 		= sys.argv[1]
my_args["id_tbl_out"] 		= sys.argv[2]
my_args["db_index_in"] 		= sys.argv[3]
my_args["db_index_out"] 	= sys.argv[4]
my_args["db_identity"] 		= sys.argv[5]

my_files = {}
my_files["name_tbl_in"] 	= "{}".format( my_args["name_tbl_in"] )
my_files["id_tbl_out"] 		= "{}".format( my_args["id_tbl_out"] )
my_files["db_index_in"] 	= "{}".format( my_args["db_index_in"] )
my_files["db_index_out"] 	= "{}".format( my_args["db_index_out"] )

my_files["db_identity_data"] 		= "{}".format( my_args["db_identity"] )
my_files["db_identity_index"] 	= "{}.index".format( my_args["db_identity"] )
my_files["db_identity_type"] 		= "{}.dbtype".format( my_args["db_identity"] )

# load name lookup table
load_name_lookup_table(my_files, my_ids)

# verify one-to-one mapping in n2i dictionare
p_keys = set(my_ids["profile"]["n2i"].keys())
q_keys = set(my_ids["query"]["n2i"].keys())
if (p_keys != q_keys):
	eprint("ERROR: Bad mapping.")
	eprint("P_KEYS:", p_keys)
	eprint("Q_KEYS:", q_keys)
	exit(0)
# verify one-to-one mapping in i2n dictionaries
p_keys = set(my_ids["profile"]["i2n"].keys())
q_keys = set(my_ids["query"]["i2n"].keys())
if (p_keys != q_keys):
	eprint("ERROR: Bad mapping.")
	eprint("P_KEYS:", p_keys)
	eprint("Q_KEYS:", q_keys)
	exit(0)

# build translation table
build_translation_table(my_files, my_ids)

# translate database index (profile -> query)
translate_db_index(my_files, my_ids)

# make identity database
make_identity_database(my_files, my_ids)

print("# completed successfully.")
