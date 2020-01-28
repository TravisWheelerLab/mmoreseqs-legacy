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
import matplotlib
import matplotlib.pyplot as plt

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

# load matrix in from .tsv file
def load_matrix(tsv_matrix):
   with open(tsv_matrix) as f:
      for line in f:
         if line[0] == "#":
            continue 

         if line.startswith("X"):
            line = line.split("\t")
            m = int(line[1])
            n = int(line[2])
            NM_MX = np.zeros((n+1,m+1,nm_st))
            SP_MX = np.zeros((m+1,sp_st))
            continue

         # all rows are tab-delimited
         line = line.split("\t")
         label = line[0].split(" ")
         l = label[0]

         # read in row from main matrix
         if l in nm_dict.keys():
            mx_cur = NM_MX
            st_cur = nm_dict[l]
            row_cur = (int)(label[1])

            for i in range(n):
               mx_cur[i][row_cur][st_cur] = (float)(line[i+1])
            continue

         # read in row from special state matrix
         if l in sp_dict.keys():
            mx_cur = SP_MX
            st_cur = sp_dict[l]

            for i in range(m+1):
               mx_cur[i][st_cur] = (float)(line[i+1])
            continue
   return NM_MX, SP_MX.T

# Find non-infinite min/max (and )
def minmax_matrix(NM_MX, SP_MX):
   min_val = np.inf
   max_val = -np.inf

   # find min and max non-infinite values
   for i in range( NM_MX.shape[0] ):
      for j in range( NM_MX.shape[1] ):
         for k in range( NM_MX.shape[2] ):
            val = NM_MX[i][j][k]
            if (val != np.inf and val != -np.inf ):
               if min_val > val:
                  min_val = val
               if max_val < val:
                  max_val = val

   # find min and max non-infinite values
   for i in range( SP_MX.shape[0] ):
      for j in range( SP_MX.shape[1] ):
         val = SP_MX[i][j]
         if (val != np.inf and val != -np.inf ):
            if min_val > val:
               min_val = val
            if max_val < val:
               max_val = val

   return min_val, max_val

# Replace non-infinite values with inf=max, -inf=min
def remove_inf_matrix(NM_MX, SP_MX, min_val, max_val):
   for i in range( NM_MX.shape[0] ):
      for j in range( NM_MX.shape[1] ):
         for k in range( NM_MX.shape[2] ):
            val = NM_MX[i][j][k]
            if (val == np.inf):
               NM_MX[i][j][k] = max_val
            elif (val == -np.inf):
               NM_MX[i][j][k] = min_val

   for i in range( SP_MX.shape[0] ):
      for j in range( SP_MX.shape[1] ):
         val = SP_MX[i][j]
         if (val == np.inf):
            SP_MX[i][j] = max_val
         elif (val == -np.inf):
            SP_MX[i][j] = min_val

# Center matrix data
def center_matrix(NM_MX, min_val, max_val, rang, abs_val):
   NM_MX -= min_val

# Scale matrix data
def scale_matrix(NM_MX, min_val, max_val, rang, abs_val):
   NM_MX *= (1 / abs_val)
   # NM_MX *= 255

# Normalize matrix
def normalize_matrix(NM_MX, min_val, max_val, rang, abs_val):
   NM_MX -= min_val
   NM_MX *= (1 / abs_val)
   # NM_MX *= 255
   # NM_MX = NM_MX.astype(int)

# Output matrix at heatmap
def output_heatmap(MAT_MX, INS_MX, DEL_MX, min_val, max_val, abs_val):
   fig, ( (ax1, ax2, ax3) ) = plt.subplots(1, 3)
   ax1.set_title( "MATCH" )
   ax1.imshow( MAT_MX, cmap='seismic', interpolation='nearest' )
   ax2.set_title( "INSERT" )
   ax2.imshow( INS_MX, cmap='seismic', interpolation='nearest' )
   ax3.set_title( "DELETE" )
   ax3.imshow( DEL_MX, cmap='seismic', interpolation='nearest' )
   # ax4.imshow( SP_MX, cmap='seismic', interpolation='nearest' )

   # center colormap
   # plt.clim(-abs_val, abs_val)

   plt.tight_layout()
   plt.show()


##############################################################################
###########################         MAIN         #############################
##############################################################################

# Import matrix
if (len(sys.argv) != 3):
   print("Usage: <tsv_matrix_1> <tsv_matrix_2>")
   sys.exit(0)
else:
   tsv_matrix = []
   tsv_matrix.append(sys.argv[1])
   tsv_matrix.append(sys.argv[2])

# number of matrices
N = 2

# Load matrix
NM_MX = []
SP_MX = []
for i in range(N):
   N_MX, S_MX = load_matrix(tsv_matrix[i])
   NM_MX.append(N_MX)
   SP_MX.append(S_MX)

# Find the difference between matrices
NM_MX = NM_MX[0] - NM_MX[1]
SP_MX = SP_MX[0] - SP_MX[1]

# Find min, max, range of matrix
min_val = np.inf
max_val = -np.inf
mx_min, mx_max = minmax_matrix(NM_MX, SP_MX)
min_val = min( min_val, mx_min )
max_val = max( max_val, mx_max )
rang = max_val - min_val

# Find the range such that color map range is centered at zero
abs_val = max( abs(min_val), abs(max_val) )


# Remove infinite values from matrix
remove_inf_matrix(NM_MX, SP_MX, min_val, max_val)


# Normalize matrix
scale_matrix(NM_MX, min_val, max_val, rang, abs_val)
scale_matrix(SP_MX, min_val, max_val, rang, abs_val)


# Split matrix into state matrices: match, delete, insert
MAT_MX = NM_MX[:,:,0] 
INS_MX = NM_MX[:,:,1] 
DEL_MX = NM_MX[:,:,2] 


# Print matrices
print( 'MAT:\n', MAT_MX )
print( 'INS:\n', INS_MX )
print( 'DEL:\n', DEL_MX )
print( 'SPECIAL:\n', SP_MX )


# Output matrices as heatmaps
output_heatmap(MAT_MX, INS_MX, DEL_MX, min_val, max_val, abs_val)

