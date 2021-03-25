 #!/usr/bin/env python
 ###############################################################################
 #     FILE:  	m8_id_appender.py
 #  PURPOSE: 	Reformats .m8 file to include 
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


# load cloud lookup file
def index_file( filename ):
	data = []

	return data


##############################################################################
###########################         MAIN         #############################
##############################################################################

# parse commandline args
pwd = os.getcwd()
if len(sys.argv) == 4:
	m8_file = sys.argv[1]
else:
   print('Usage: <index_file>')
   sys.exit(0)


# m8 file
with open( m8_file, "r" ) as fp:
	for readline in fp:
		pass