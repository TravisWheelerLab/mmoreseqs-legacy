#!/usr/bin/env python

###############################################################################
#  @file prob_to_bitmap
#  @brief Takes in tsv of dp matrix and creates a bitmap.
#
#  @author Dave Rich
#  @bug Lots.
###############################################################################

import sys
import os
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
   edg["FILE"] = filename;
   edg["DATA"] = {}
   fp = open( filename, "r" )

   for line in fp:
      line = line.replace(","," ")
      
      if line.startswith("#"):
         fields = line.split()
         if fields[1].startswith("N:"):
            edg["N"] = fields[2]
         if fields[1].startswith("ORIENTATION:"):
            edg["ORIENTATION"] = fields[2]
         if fields[1].startswith("Q="):
            Q = int(fields[1].split("=")[1])
            T = int(fields[2].split("=")[1])
            edg["Q"] = Q
            edg["T"] = T
      
      if line.startswith("["):
         fields   = line.split()
         row      = int(fields[3])
         my_lb    = int(fields[5])
         my_rb    = int(fields[7])

         edg["DATA"][row] = ( my_lb, my_rb )

   return edg

# Add edgebounds to matrix
def add_edgebounds_to_matrix( mx, edg ):
   for row in edg["DATA"].keys():
      bounds = edg["DATA"][row]
      lb = bounds[0]
      rb = bounds[1]
      
      for col in range(lb, rb):
         mx[row,col] += 1
   return

# Output heatmap of single matrix with traceback on top
def output_heatmap( mx ):
   fig, ax1 = plt.subplots(1, 1)

   # plot heatmap
   # ax1.set_title( title )
   im = ax1.imshow( mx, cmap='jet' )

   # # plot viterbi traceback
   # tr_len = len(TR[0])-1

   # # draw the viterbi trace
   # if tr_len > 2:
   #      ax1.scatter( TR[0][0], TR[1][0], c='white', s=3 )
   #      ax1.scatter( TR[0][tr_len], TR[1][tr_len], c='white', s=3 )

   #      #ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')
   #      ax1.plot( TR[0], TR[1], linestyle='-', linewidth=2, color='white')
   #      ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')

   if save_fig:
      dest = "{}/{}.LAYMAP.jpg".format(out_file, name)
      plt.savefig(dest)
      print("Figure saved to '{}'...".format(dest))

   if show_fig:
      plt.show()
   
   return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# default location to save files
out_file = os.getcwd()
in_files = []
name = "layer-map"
show_fig = True
save_fig = True

# Parse args
if (len(sys.argv) <= 1):
   print("Usage: <edg_file_1> <...>");
   sys.exit(1);
else:
   i = 1
   while (i < len(sys.argv) ):
      arg = sys.argv[i]
      if arg.startswith("--"):
         i += 1
         opt = sys.argv[i]
         if arg == "--name":
            name = opt 
         if arg == "--no-show":
            show_fig = False
      else:
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
mx = build_matrix(Q+1,T+1)

# add edgebounds 
for i in range(len(edges)):
   edg = edges[i]
   add_edgebounds_to_matrix( mx, edg )

# output results
output_heatmap( mx )
