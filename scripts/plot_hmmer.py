 #!/usr/bin/env python
 ###############################################################################
 #  @file  		filter_results.py
 #  @brief 		Filter results
 #
 #  @author 	Dave Rich
 #  @bug 		Lots.
 ###############################################################################

import sys
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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


# load lookup
def load_lookup( filename ):
	data = {}
	return data


##############################################################################
###########################         MAIN         #############################
##############################################################################

# get args
if ( len(sys.argv) <= 1 ):
	print("Usage: <hmmer_times>")
	exit(0)

hmmer_time_fname = sys.argv[1]

# load data
hmmer_data = load_hmmer_time(hmmer_time_fname)

# linear spacing
x = range( len(mmseqs_plus_data) )

# sort data by size
hmmer_data.sort( key=lambda r:int(r[7]) )

# plot time data against problem size
fwdbck_time = []
size = []
for data in hmmer_data:
	fwdbck_time.append( float(data[4]) )
	size.append( int(data[7]) )

plt.scatter( size, fwdbck_time, marker='o', s=2 )
plt.show()