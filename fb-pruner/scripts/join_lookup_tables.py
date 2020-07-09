 #!/usr/bin/env python
 ###############################################################################
 #     FILE:  	sync_lookup_ids.py
 #  PURPOSE: 	Pairs mmseqs lookup ID with cloudsearch lookup ID.
 # 				Outputs as ID pairs.
 #
 #   AUTHOR: 	Dave Rich
 #     BUGS: 	
 ###############################################################################

import sys, os
import numpy as np
# import cv2 as cv
# from PIL import Image
# import matplotlib
# import matplotlib.pyplot as plt

# load mmseqs lookup file
def load_mmseqs_lookup( filename, ftype = "dict" ):
	data = {}

	with open( filename, "r" ) as fp:
		for line in fp:
			# skip comment lines
			if line.startswith("#"):
				continue
			# check for proper number of fields
			fields = line.split()
			if len(fields) < 3:
				continue
			# get data fields
			mid  = fields[0]
			name = fields[1]
			offset = fields[2]

			#insert into database
			data[name] = mid

	return data

# load cloud lookup file
def load_cloud_lookup( filename, ftype = "dict" ):
	data = {}

	with open( filename, "r" ) as fp:
		for line in fp:
			# skip comment lines
			if line.startswith("#"):
				continue
			# check for proper number of fields
			fields  = line.split()
			if len(fields) < 3:
				continue
			# get data fields
			cid	    = fields[0]
			offset 	= fields[1]
			name 	= fields[2]

			# insert into database
			data[name] = cid

	return data

# join mmseqs and cloud lookup tables
def join_lookup_dicts( mmseqs_dict, cloud_dict ):
	joint_dict = {}
	mmseqs_names = set( mmseqs_dict.keys().sort() )
	cloud_names = set( cloud_dict.keys().sort() )
	for name in mmseqs_names:
		if name in cloud_names:
			joint_dict[name] = ( cloud_dict[name], mmseqs_dict[name] )
		else:
			print("# ERROR: '{}' from mmseqs_dict not found in cloud_dict.".format(name) )

	return joint_dict

# load joint lookup file
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
			data.append( (cid, offset, name) )

	# sort ids according to names 
	data.sort( key = lambda x: x[2] )
	return data

# check if names have the same root
def rootname_match( mmseqs_name, cloud_name ):
	min_len = min( len(cloud_name), len(mmseqs_name) )
	same_name = ( cloud_name[0:min_len] == mmseqs_name[0:min_len] )
	return same_name

def print_joint_data( joint_data ):
	# print header
	print('#cloud_id\t #mmseqs_id\t #cloud_name')

	# pair ids from each 
	for name in joint_data.keys():

		data = joint_data[name]
		print( "{}\t{}\t{}".format( data[0], data[1], name ) )

	return

##############################################################################
###########################         MAIN         #############################
##############################################################################

# parse commandline args
pwd = os.getcwd()
if len(sys.argv) == 3:
	lookup_mmseqs = sys.argv[1]
	lookup_cloud  = sys.argv[2]	
else:
   print('Usage: <lookup_mmseqs> <lookup_cloud> <output_file>')
   sys.exit(0)

# load data
print("# load mmseqs...")
mmseqs_data = load_mmseqs_lookup( lookup_mmseqs )
print("# load cloud....")
cloud_data = load_cloud_lookup( lookup_cloud )
print("# join data...")
joint_data = join_lookup_dicts( mmseqs_data, cloud_data )

# number of proper matches
match_cnt = 0
# get length 
N = len(mmseqs_data)

print_joint_data( joint_data )