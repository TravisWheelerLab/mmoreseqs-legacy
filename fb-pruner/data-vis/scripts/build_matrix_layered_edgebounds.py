#!/usr/bin/env python

###############################################################################
#  @file prob_to_bitmap
#  @brief Takes in tsv of dp matrix and creates a bitmap.
#
#  @author Dave Rich
#  @bug Lots.
###############################################################################

import sys
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# set resolution of output
mpl.rcParams['figure.dpi'] = 200

# constants: normal states, special states
nm_st = 3
sp_st = 5

nm_dict = {
   "M": 0, 
   "I": 1, 
   "D": 2
}

sp_dict = {
   "N": 0, 
   "J": 1, 
   "E": 2, 
   "C": 3, 
   "B": 4
}

tsv_matrix = [];

# build matrix
def build_matrix( m, n ):
   mx = np.zeros((m,n))
   return mx

# Load edgebounds from file
def load_edg( filename ):
   edg = {}
   fp = open( filename, "r" )

   for line in fp:
      line = line.replace(","," ")
      
      if line.startswith("#"):
         fields = line.split()
         if fields[1].startswith("N:"):
            edg["N"] = fields[2]
         if fields[1].startswith("ORIENTATION:"):
            edg["orientation"] = fields[2]
         if fields[1].startswith("Q="):
            Q = int(fields[1].split("=")[1])
            T = int(fields[2].split("=")[1])
            edg["Q"] = Q
            edg["T"] = T
      
      if line.startswith("["):
         fields = line.split()
         my_id = int(fields[3])
         my_lb = int(fields[5])
         my_rb = int(fields[5])

         edg[my_id] = ( my_lb, my_rb )

   print(edg)

   return edg

# Add edgebounds to matrix
def add_edgebounds_to_matrix( mx, edges ):
   
   for my_id in edges.keys():
      bounds = edges[my_id]
      lb = bounds[0]
      rb = bounds[1]
      for i in range(lb, rb):
         mx[my_id, i] += 1

   return

# Output heatmap of single matrix with traceback on top
def output_heatmap( mx ):
   fig, ax1 = plt.subplots(1, 1)

   # plot heatmap
   # ax1.set_title( title )
   im = ax1.imshow( mx, cmap='jet', interpolation='nearest' )

   # # plot viterbi traceback
   # tr_len = len(TR[0])-1

   # # draw the viterbi trace
   # if tr_len > 2:
   #      ax1.scatter( TR[0][0], TR[1][0], c='white', s=3 )
   #      ax1.scatter( TR[0][tr_len], TR[1][tr_len], c='white', s=3 )

   #      #ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')
   #      ax1.plot( TR[0], TR[1], linestyle='-', linewidth=2, color='white')
   #      ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')

   # if (opts["-S"]):
   #    dest = "{}/{}.TR_FIG.jpg".format(file, name)
   #    plt.savefig(dest)
   #    print("Figure saved to '{}'...".format(dest))

   plt.show()
   return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# default location to save files
out_file = "/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data-vis/heatmaps/"
in_files = []

# Parse args
if (len(sys.argv) <= 1):
   print("Usage: <edg_file_1> <...>");
   sys.exit(1);
else:
   i = 1
   while (i < len(sys.argv) ):
      in_files.append(sys.argv[i])
      i += 1

edges = []

# load in edgebounds
for i in range(len(in_files)):
   print("{}: {}".format(i, in_files[i]))
   in_file = in_files[i]
   edg = load_edg(in_file)
   edges.append(edg)
   pass

# verify all edgebounds
Q = edges[i]["Q"]
T = edges[i]["T"] 
for i in range(len(edges)):
   edg = edges[i]
   if (Q != edg["Q"]) or (T != edg["T"]):
      print("Error: not all edges are the same dimension.")
      sys.exit(1)
   pass

# create matrix
mx = build_matrix(Q,T)

# add edgebounds 
for i in range(len(edges)):
   edg = edges[i]
   add_edgebounds_to_matrix( mx, edg )

# output results
output_heatmap( mx )
