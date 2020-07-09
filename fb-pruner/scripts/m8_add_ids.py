 #!/usr/bin/env python
 ###############################################################################
 #     FILE:  	m8_id_appender.py
 #  PURPOSE: 	Reformats .m8 file to include ids.
 # 				Adds target id and query id fields to .m8 file and prints results. 
 #
 #   AUTHOR: 	Dave Rich
 #     BUGS: 	
 ###############################################################################

import sys, os
import numpy as np

# load mmseqs lookup from file
def load_mmseqs_lookup( filename ):
	data = {}

	with open( filename, "r" ) as fp:
		for line in fp:
			# skip comment lines
			if line.startswith("#"):
				continue
			# check for proper number of fields
			line = line.split()
			if len(line) < 3:
				continue
			# get data fields
			mid  = line[0]
			name = line[1]
			#insert into list
			data[name] = mid

	return data

# load cloud lookup from file
def load_cloud_lookup( filename ):
	data = []

	with open( filename, "r" ) as fp:
		for line in fp:
			# skip comment lines
			if line.startswith("#"):
				continue
			# check for proper number of fields
			line 	= line.split()
			if len(line) < 3:
				continue
			# get data fields
			tid	    = line[0]
			offset 	= line[1]
			name 	= line[2]
			# insert into dictionary
			data.append( (tid, offset, name) )

	# sort ids according to names 
	data.sort( key = lambda x: x[2] )
	return data

# load joint lookup from file
def load_joint_lookup( filename ):
	data = []

	with open( filename, "r" ) as fp:
		for line in fp:
			# skip comment lines
			if line.startswith("#"):
				continue
			# check for proper number of fields
			line 	= line.split()
			if len(line) < 3:
				continue
			# get data fields
			cid	    = line[0]
			mid 	= line[1]
			name 	= line[2]
			# insert into dictionary
			data.append( (cid, mid, name) )

	# sort ids according to names 
	# data.sort( key = lambda x: x[2] )
	return data

def add_ids_to_m8( m8_file, t_lookup, q_lookup ):
	res_id = 0

	# m8 file
	with open( m8_file, "r" ) as fp:
		for line in fp:
			fields = line.split()

			# extract fields from result
			q_name 	= fields[0]
			t_name 	= fields[1]

			# look up the ids
			t_id = t_lookup[t_name]
			q_id = q_lookup[q_name]

			line = line.strip()

			# append result id, query id, target id to front of results line
			print( "{}\t{}\t{}\t{}\t{}\t{}".format(res_id, q_id, t_id, q_id, t_id, line) )

			res_id += 1

	return



##############################################################################
###########################         MAIN         #############################
##############################################################################

# parse commandline args
pwd = os.getcwd()
if len(sys.argv) == 4:
	m8_file = sys.argv[1]
	t_lookup_file = sys.argv[2]
	q_lookup_file = sys.argv[3]
else:
   print('Usage: <m8_results_file> <t_lookup> <q_lookup>')
   sys.exit(0)

t_lookup = load_mmseqs_lookup(t_lookup_file)
q_lookup = load_mmseqs_lookup(q_lookup_file)

add_ids_to_m8( m8_file, t_lookup, q_lookup )
