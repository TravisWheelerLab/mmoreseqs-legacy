 #!/usr/bin/env python
 ###############################################################################
 #  @file  analyze_results
 #  @brief Render benchmark scores as a line graph
 #
 #  @author Dave Rich
 #  @bug Lots.
 ###############################################################################

import sys
import numpy as np
import pandas as pd
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# loads .mym8 file, fwdbck is True if fwdbck
def load_mym8( filename ):
 	# field counter
 	f_cnt = 0
	data = []

 	with open( filename, "r" ) as fp:
 		for line in fp:
 			# check for comment lines or empty lines
 			if ( (line[0] == ">") or (line[0] == "#") ):
 				continue

 			# break up line
 			fields = line.split()

 		# 	# ids
 		# 	result_id 	= fields[f_cnt]
 		# 	t_id 		= fields[1]
 		# 	q_id  		= fields[2]

 		# 	# problem size
 		# 	t_len		= fields[3]
 		# 	q_len		= fields[4]
 		# 	total_cells = fields[5]
 		# 	cloud_cells = fields[6]

 		# 	# parameters
 		# 	alpha 		= fields[7]
 		# 	beta 		= fields[8]

 		# 	# scores
		# 	vit_sc 		= fields[9]
		# 	fwd_sc 		= fields[10]
		# 	bck_sc 		= fields[11]
 		# 	cld_sc 		= fields[12]

 		# 	# times
		# 	vit_t 		= fields[13]
		# 	fwd_t 		= fields[14]
		# 	bck_t 		= fields[15]

		# 	# cloud time
 		# 	cld_fwd_t 	= fields[16]
 		# 	cld_bck_t	= fields[17]
 		# 	merge_t 	= fields[18]
 		# 	reorient_t 	= fields[19]
 		# 	bnd_fwd_t 	= fields[20]
 		# 	total_t 	= fields[21]

 			data.append(fields)
 			continue

 	return data

# load hmmer 
def load_hmmer_time( filename ):
	# field counter
	f_cnt = 0
	data = []

	with open( filename, "r" ) as fp:
		for line in fp:
			# check for comment lines or empty lines
 			if ( line.startswith("##") ):
 				continue
 			fields 		= line.split()
 			# check that there are a proper number of fields 
 			if ( len(fields) != 7 ):
 				continue

 			total_cells = ( int(fields[0])+1 ) * ( int(fields[2])+1 )
 			fields.append( total_cells )

 			data.append( fields )
 			continue
 	return data

# # get the length of each bin on the y-axis data
# def bin_lengths( data_y, bin_list ):
# 	bin_lengths = []
# 	# track length of current bin
# 	cnt = 0
# 	# track bin list element
# 	j = 0 
# 	for i in range( len(data_y) ):
# 		if (data[i] > bin_list[j] )


# # bin data on the x-axis using bin length of y-axis
# def bin_data( data_x, bin_length ):
# 	binned_data = []
# 	i = 0
# 	for val in bin_list:
# 		tot = 0
# 		cnt = 0
# 		while( data[i] < bin_list ):
# 			tot += data[i]
# 			cnt += 1



##############################################################################
###########################         MAIN         #############################
##############################################################################

if ( len(sys.argv) <= 1 ):
	print( "Usage: <mmseqs_plus_fname> <hmmer_time_fname>" )
	exit(0)

# load data
mmseqs_plus_fname  	= sys.argv[1]
mmseqs_plus_data = load_mym8(mmseqs_plus_fname)

# linear spacing
x = range( len(mmseqs_plus_data) )

# load data
# hmmer_time_fname	= sys.argv[2]
# hmmer_data = load_hmmer_time(hmmer_time_fname)

# sort data by total_cells (problem size)
mmseqs_plus_data.sort( key=lambda r:int(r[5]) )

# capture problem size data
cloud_size = []
size = []
for data in mmseqs_plus_data:
	size.append( int(data[5]) )
	cloud_size.append( int(data[6]) )

# capture time data
cloud_time = []
fwdbck_time = []
fwd_time = []
vit_time = []
for data in mmseqs_plus_data:
	vit_time.append( float(data[13]) )
	fwd_time.append( float(data[14]) )
	fwdbck_time.append( float(data[14]) + float(data[15]) )
	cloud_time.append( float(data[21]) )
times = []

# 
for i in range(len(vit_time)):
	times.append( [ fwdbck_time[i], fwd_time[i], vit_time[i], cloud_time[i] ] )

# plt.scatter( size, cloud_time, marker='o', s=2, color='red' )
plt.scatter( size, fwdbck_time, marker='o', s=2, color='red' )
plt.show()


# plot score data against problem size
# cloud_score = []
# fwdbck_score = []
# fwd_score = []
# vit_score = []
# scores = [ vit_score, fwd_score, fwdbck_score, cloud_score ]
# for i in len(data):
# 	vit_time.append( mmseqs_plus_data[i][9] )
# 	fwd_score.append( mmseqs_plus_data[i][10] )
# 	fwdbck_score.append( mmseqs_plus_data[i][10] + mmseqs_plus_data[i][11] )
# 	cloud_score.append( mmseqs_plus_data[i][12] )
# plt.plot( times, size, labels=['Cloud','Forward','Forward-Backward'])
# plt.legend(loc='upper left')
# plt.show()