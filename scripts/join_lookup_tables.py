 #!/usr/bin/env python
 ###############################################################################
 #     FILE:  	sync_lookup_ids.py
 #  PURPOSE: 	Pairs mmseqs lookup ID with cloudsearch lookup ID.
 # 				Outputs as ID pairs
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

# parse commandline args
pwd = os.getcwd()
if len(sys.argv) == 3:
	lookup_mmseqs = sys.argv[1]
	lookup_cloud  = sys.argv[2]	
else:
   print('Usage: <lookup_mmseqs> <lookup_cloud>')
   sys.exit(0)

# load mmseqs lookup file
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
			mid  = line[0]
			name = line[1]
			#insert into list
			data.append( (mid, name) )

	# sort ids according to names 
	data.sort( key = lambda x: x[1] )
	return data

# load cloud lookup file
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
			cid	    = line[0]
			offset 	= line[1]
			name 	= line[2]
			# insert into dictionary
			data.append( (cid, offset, name) )

	# sort ids according to names 
	data.sort( key = lambda x: x[2] )
	return data

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

##############################################################################
###########################         MAIN         #############################
##############################################################################

# load data
mmseqs_data = load_mmseqs_lookup( lookup_mmseqs )
cloud_data = load_cloud_lookup( lookup_cloud )
# print( "mmseqs_data_len: ", len(mmseqs_data) )
# print( "cloud_data_len:", len(cloud_data) )

# number of proper matches
match_cnt = 0
# get length 
N = len(mmseqs_data)

# print header
# print('#cloud_id\t #mmseqs_id\t #cloud_name')

# pair ids from each 
for i in range(N):
	cloud  	 = cloud_data[i]
	mmseqs 	 = mmseqs_data[i]
	cloud_id, cloud_offset, cloud_name = cloud
	mmseqs_id, mmseqs_name = mmseqs

	# test check if mmseqs and cloud names have the same root
	# name_check = rootname_match( mmseqs_name, cloud_name )
	# if same_name == False:
	# 	print(f'{i}\t {cloud_name}\t {mmseqs_name}\t')

	# format; cloud_id, cloud_offset, (cloud_name), mmseqs_id, mmseqs_name
	# print(f'{cloud_id}\t {cloud_offset}\t {cloud_name}\t {mmseqs_id}\t {mmseqs_name}')
	# print(f'{cloud_name}\t {mmseqs_name}\t')
	print(f'{cloud_id}\t {mmseqs_id}\t {mmseqs_name}')

	continue
