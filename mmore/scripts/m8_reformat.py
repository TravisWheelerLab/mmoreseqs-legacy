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
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt


# load mmseqs lookup from file
def load_mmseqs_lookup( filename ):
	data = []

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
			tid  = line[0]
			name = line[1]
			#insert into list
			data.append( (tid, name) )

	# sort ids according to names 
	data.sort( key = lambda x: x[1] )
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

# create lookup based
def create_joint_lookup_dict( joint_data ):
	name_to_id = {}

	for i in range( len(joint_data) ):
		cid, mid, name = joint_data[i]
		name_to_id[ name ] = ( cid, mid ) 
		pass

	return name_to_id


##############################################################################
###########################         MAIN         #############################
##############################################################################

# parse commandline args
pwd = os.getcwd()
if len(sys.argv) == 4:
	m8_file = sys.argv[1]
	q_joint_lookup = sys.argv[2]
	t_joint_lookup = sys.argv[3]
else:
   print('Usage: <m8_results_file> <q_joint_lookup> <t_joint_lookup>')
   sys.exit(0)

# print header
# print(f'#        m8_file:\t{m8_file}')
# print(f'# q_joint_lookup:\t{q_joint_lookup}')
# print(f'# t_joint_lookup:\t{t_joint_lookup}')

# load data
q_joint_data = load_joint_lookup( q_joint_lookup )
t_joint_data = load_joint_lookup( t_joint_lookup )
q_name_to_id = create_joint_lookup_dict( q_joint_data )
t_name_to_id = create_joint_lookup_dict( t_joint_data )

res_id = 0

# m8 file
with open( m8_file, "r" ) as fp:
	for readline in fp:
		line = readline.split()

		# extract fields from result
		q_name 	= line[0]
		t_name 	= line[1]
		# perc_id 	= line[2]
		# aln_len 	= line[3]
		# mismatch  = line[4]
		# gap_open  = line[5]
		# q_start 	= line[6]
		# q_end 	= line[7]
		# t_start 	= line[8]
		# t_end 	= line[9]
		# e_val 	= line[10]
		# bit_sc 	= line[11]

		# look up the ids
		t_cid, t_mid = t_name_to_id[ t_name ]
		q_cid, q_mid = q_name_to_id[ q_name ]

		# append result id, query id, target id to front of results line
		print(f'{res_id}\t', end="")
		print(f'{q_cid}\t', end="")
		print(f'{t_cid}\t', end="")
		print(f'{q_mid}\t', end="")
		print(f'{t_mid}\t', end="")
		print(readline, end="")

		res_id += 1
