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
   edges = {}
   fp = open( filename, "r" )

   for line in fp:
      if line.startswith("#"):
         fields = line.split()
      else:
         fields = line.split()
         my_id = int(fields[3].replace(",",""))
         my_lb = int(fields[5].replace(",",""))
         my_rb = int(fields[5].replace(",",""))

         edg[my_id] = [ my_lb, my_rb ]

   return edges

# Add edgebounds to matrix
def add_edgebounds_to_matrix( mx, edges ):
   
   for my_id in edges.keys():
      bounds = edges[my_id]

      for i in range(lb, rb):
         mx[my_id, ]

   return

# Output matrix at heatmap
def output_heatmap(title, MAT_MX, INS_MX, DEL_MX, SP_MX, vmin, vmax, file="", default="test"):
   fig, ( (ax1, ax2), (ax3, ax4) ) = plt.subplots(2, 2)

   ax1.set_title( title )
   ax1.imshow( MAT_MX, cmap='jet', interpolation='nearest' )
   ax2.set_title( "INSERT" )
   ax2.imshow( INS_MX, cmap='jet', interpolation='nearest' )
   ax3.set_title( "DELETE" )
   ax3.imshow( DEL_MX, cmap='jet', interpolation='nearest' )
   ax4.set_title( "SPECIAL" )
   im = ax4.imshow( SP_MX, cmap='jet', interpolation='nearest' )

   # Add legend
   cbar = ax1.figure.colorbar(im, ax=ax4)

   plt.tight_layout()
   if (opts["-S"]):
      dest = "{}/{}.TR_FIG.jpg".format(file, name)
      plt.savefig(dest)
      print("Figure saved to '{}'...".format(dest))

   plt.show()
   return None

# Output heatmap of single matrix with traceback on top
def output_heatmap_trace(title, MAT_MX, TR, vmin, vmax, file="", name="test"):
   fig, ax1 = plt.subplots(1, 1)

   # plot heatmap
   ax1.set_title( title )
   im = ax1.imshow( MAT_MX, cmap='jet', interpolation='nearest' )

   # plot viterbi traceback
   tr_len = len(TR[0])-1

   # draw the viterbi trace
   if tr_len > 2:
        ax1.scatter( TR[0][0], TR[1][0], c='white', s=3 )
        ax1.scatter( TR[0][tr_len], TR[1][tr_len], c='white', s=3 )

        #ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')
        ax1.plot( TR[0], TR[1], linestyle='-', linewidth=2, color='white')
        ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')

      # Create a Rectangle patch for Viterbi window
      # rect = patches.Rectangle((50,100),40,30,linewidth=1,edgecolor='r',facecolor=None)
      # ax1.add_patch(rect)

      # Add legend
      #cbar = ax1.figure.colorbar(im, ax=ax1)

   if (opts["-S"]):
      dest = "{}/{}.TR_FIG.jpg".format(file, name)
      plt.savefig(dest)
      print("Figure saved to '{}'...".format(dest))

   plt.show()
   return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# default location to save files
file = "/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data-vis/heatmaps/"

# Parse args
if (len(sys.argv) == 1):
   print("Usage: <tsv_matrix_1>");
   sys.exit(0);
else:
   i = 1
   while (i < len(sys.argv) ):
      arg = sys.argv[i];
      # apply options
      if (arg.startswith("-")):
         if arg in opts.keys():
            opts[arg] = True;
      else:
         tsv_matrix.append(sys.argv[i]);
      i += 1

