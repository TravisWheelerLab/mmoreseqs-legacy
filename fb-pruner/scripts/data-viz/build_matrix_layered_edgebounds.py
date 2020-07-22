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
import re

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

         if row not in edg["DATA"]:
            edg["DATA"][row] = []
         edg["DATA"][row].append( ( my_lb, my_rb ) )

   return edg

# Add edgebounds to matrix
def add_edgebounds_to_matrix( mx, edg ):
   for row in edg["DATA"].keys():
      bounds = edg["DATA"][row]

      for bound in bounds:
         lb = bound[0]
         rb = bound[1]
      
         for col in range(lb, rb):
            mx[row,col] += 1
   return

# count cells
def cell_count( edg ):
   count = 0
   for row in edg["DATA"].keys():
      bounds = edg["DATA"][row]
      lb = bounds[0]
      rb = bounds[1]
      count += ( rb - lb )
   return count

# Load trace file
def load_trace( trace_file ):
   TR = [[],[]]

   with open( trace_file, "r" ) as fp:
      for line in fp:
         if line.startswith("#"):
            continue
         fields = line.replace("(", ",").replace(")", ",").split(",")
         TR[0].append( fields[2] )
         TR[1].append( fields[3] )

   return TR

# Output heatmap of single matrix with traceback on top
def output_heatmap( mx ):
   mx = mx.T
   fig, ax1 = plt.subplots(1, 1)

   # plot heatmap
   # ax1.set_title( title )
   im = ax1.imshow( mx, cmap='jet' )

   # plot viterbi traceback
   if trace != None:
      tr_end = len(trace[0])-1

      # draw the viterbi trace
      if tr_len > 2:
           ax1.scatter( trace[0][0], trace[1][0], c='white', s=3 )
           ax1.scatter( trace[0][tr_end], trace[1][tr_end], c='white', s=3 )

           ax1.plot( trace[0], trace[1], linestyle='-', linewidth=2, color='white')
           ax1.plot( trace[0], trace[1], linestyle='-', linewidth=1, color='black')


   plt.xlabel( t_name )
   plt.ylabel( q_name )
   plt.title( title )

   if save_fig:
      dest = "{}/{}".format(out_folder, out_file)
      plt.savefig(dest)
      print("Figure saved to '{}'...".format(dest))

   if show_fig:
      plt.show()
   
   return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# default location to save files
out_folder = os.getcwd()
in_files = []
out_file = "layer-map.LAYMAP.jpg"

show_fig = True
save_fig = True
trace_file = None
trace = None

title = "Alpha/Beta Heatmap"
t_name = "TARGET"
q_name = "QUERY"

# Parse args
if (len(sys.argv) <= 1):
   print("Usage: <edg_file_1> <...>");
   sys.exit(1);
else:
   i = 1
   while (i < len(sys.argv) ):
      arg = sys.argv[i]
      if arg.startswith("--"):
         opt = sys.argv[i]
         if arg == "--output":
            i += 1
            out_file = sys.argv[i]
         if arg == "--title":
            i += 1
            title = sys.argv[i]
         if arg == "--trace":
            i += 1
            trace_file = sys.argv[i] 
         if arg == "--tname":
            i += 1
            t_name = sys.argv[i]
         if arg == "--qname":
            i += 1
            q_name = sys.argv[i]
         if arg == "--no-show":
            show_fig = False
         if arg == "--no-save":
            save_fig = False

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

# load in trace
if trace_file != None:
   trace = load_trace( trace_file )

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
