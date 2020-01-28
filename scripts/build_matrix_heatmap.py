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
import matplotlib.patches as patches

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

opts = {
   "-4": False,
   "-tr": False,
   "-diff": False,
   "-add": False,
   "-eq": False
};

tsv_matrix = [];


# If value is infinity, replace it with another value
def replace_inf(val, replace):
   if val == np.inf:
      return replace
   if val == -np.inf:
      return -replace
   return val

# Load matrix in from .tsv file
def load_matrix(tsv_matrix):
   trace = [[],[]]

   # read every line in file
   with open(tsv_matrix) as f:
      line = f.readline()
      while line:

         # ignore comment lines
         if line[0] == "#":
            pass

         # dimensions of matrix (single line)
         elif line.startswith("XDIM"):
            line = line.split("\t")
            m = int(line[1])
            n = int(line[2])
            NM_MX = np.zeros((n+1,m+1,nm_st))
            SP_MX = np.zeros((m+1,sp_st))
            pass

         # traceback indices
         elif line.startswith("XTRACE"):
            line = f.readline()
            while line and not line.startswith("/"):
               line = line.split("\t")
               trace[0].append( int(line[1]) )
               trace[1].append( int(line[2]) )

               line = f.readline()
            pass

         # matrix values
         elif line.startswith("XMATRIX"):
            line = f.readline()
            while line and not line.startswith("/"):
               # label which matrix it goes to
               line = line.split("\t")
               label = line[0].split(" ")
               l = label[0]

               # read in row from main matrix
               if l in nm_dict.keys():
                  mx_cur = NM_MX
                  st_cur = nm_dict[l]
                  row_cur = (int)(label[1])

                  for i in range(n+1):
                     mx_cur[i][row_cur][st_cur] = (float)(line[i+1])
                  pass

               line = f.readline()
            pass

         elif line.startswith("XSPECIAL"):
            line = f.readline()
            while line and not line.startswith("/"):
               # label which row 
               line = line.split("\t")
               label = line[0]

               # read in row from special state matrix
               if label in sp_dict.keys():
                  mx_cur = SP_MX
                  st_cur = sp_dict[label]

                  for i in range(m+1):
                     mx_cur[i][st_cur] = (float)(line[i+1])
                  pass

               line = f.readline()
            pass

         line = f.readline()

   return NM_MX, SP_MX.T, trace

# Take cell-wise difference of two matrices
def diff_matrix(mat1, mat2): 
   res = np.zeros_like(mat1)

   for i,x in np.ndenumerate(mat1):
      mat1[i] = replace_inf(mat1[i], 50)
      mat2[i] = replace_inf(mat2[i], 50)

      if (mat1[i] == mat2[i]):
         res[i] = 0
      else:
         res[i] = mat1[i] - mat2[i]

   return res

# Take cell-wise difference of two matrices
def add_matrix(mat1, mat2):
   res = np.zeros_like(mat1)

   for i,x in np.ndenumerate(mat1):
      if (mat1[i] == -mat2[i]):
         res[i] = 0
      else:
         res[i] = mat1[i] + mat2[i]

   return res

# Test equality of two matrices (within tolerance)
def eq_matrix(mat1, mat2, tol = 0.01):
   res = np.zeros_like(mat1)

   for i,x in np.ndenumerate(mat1):
      if abs(mat1[i] - mat2[i]) < tol :
         res[i] = 0
      else:
         res[i] = 1

   return res

# Find non-infinite min/max
def minmax_matrix(NM_MX, SP_MX):
   min_val = np.inf
   max_val = -np.inf

   # find min and max non-infinite values
   for i,x in np.ndenumerate(NM_MX):
      val = NM_MX[i]
      if (val != np.inf and val != -np.inf ):
         if min_val > val:
            min_val = val
         if max_val < val:
            max_val = val

   # find min and max non-infinite values
   for i,x in np.ndenumerate(SP_MX):
      val = SP_MX[i]
      if (val != np.inf and val != -np.inf ):
         if min_val > val:
            min_val = val
         if max_val < val:
            max_val = val

   return min_val, max_val

# Replace non-infinite values with inf=max*mult, -inf=min*mult
def remove_inf_matrix(NM_MX, SP_MX, min_val, max_val):
   multiplier = 1.2
   # remove infinite values
   for i,x in np.ndenumerate(NM_MX):
      val = NM_MX[i]
      if (val == np.inf):
         NM_MX[i] = max_val * multiplier
      elif (val == -np.inf):
         NM_MX[i] = min_val * multiplier

   # remove infinite values
   for i,x in np.ndenumerate(SP_MX):
      val = SP_MX[i]
      if (val == np.inf):
         SP_MX[i] = max_val * multiplier
      elif (val == -np.inf):
         SP_MX[i] = min_val * multiplier
   return None

# Normalize matrix
def normalize_matrix(NM_MX, min_val, max_val, rang):
   NM_MX -= min_val
   NM_MX *= (1 / rang)
   NM_MX = np.exp(NM_MX)
   # NM_MX *= 255
   # NM_MX = NM_MX.astype(int)
   return None

# Normalize all anti-diagonals of matrix
def normalize_antidiags(NM_MX):
   print(NM_MX.shape)
   x,y,C = NM_MX.shape
   min_corner = min(x,y)
   max_corner = max(x,y)
   num_diags = x+y-1

   # for each state (M,I,D)
   for c in range(C):
      num_cells = 0
      start_i = 0
      M = NM_MX

      print("num_diags:",num_diags)

      # for each diagonal in matrix...
      for d in range(num_diags):
         if d < max_corner:
            num_cells += 1
         if d > min_corner-1:
            num_cells -= 1
         if d > y-1:
            start_i += 1
         cell_cnt = 0
         min_val = np.inf
         max_val = -np.inf

         # find min and max values
         for i in range(start_i, start_i+num_cells):
            for c in range(C):
               j = d-i
               # M[i,j,c] = d
               cell_cnt += 1
               # print(d,i,j)
               if M[i,j,c] == np.inf or M[i,j,c] == -np.inf:
                  next
               if M[i,j,c] > max_val:
                  max_val = M[i,j,c]
               if M[i,j,c] < min_val:
                  min_val = M[i,j,c]

         # replace inf vals with min/max vals
         for i in range(start_i, start_i+num_cells):
            for c in range(C):
               j = d-i
               # M[i,j,c] = d
               cell_cnt += 1
               # print(d,i,j)
               if M[i,j,c] == np.inf:
                  M[i,j,c] = max_val
               if M[i,j,c] == -np.inf:
                  M[i,j,c] = min_val

         range_val = max_val - min_val
         # print("PRE d:",d,"range_val:",range_val,"min_val:",min_val, "max_val:", max_val)

         # normalize the diagonal
         for i in range(start_i, start_i+num_cells):
            for c in range(C):
               j = d-i
               if (range_val != 0):
                  M[i,j,c] = ((M[i,j,c] - min_val)/ range_val)
               else:
                  M[i,j,c] = 0.0
         # print("POST d:",d,"range_val:",range_val,"min_val:",min_val, "max_val:", max_val)
   return None       

# Output matrix at heatmap
def output_heatmap(title, MAT_MX, INS_MX, DEL_MX, SP_MX, vmin, vmax):
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
   plt.show()
   return None

# Output heatmap of single matrix with traceback on top
def output_heatmap_trace(title, MAT_MX, TR, vmin, vmax):
   fig, ax1 = plt.subplots(1, 1)

   # plot heatmap
   ax1.set_title( title )
   im = ax1.imshow( MAT_MX, cmap='jet', interpolation='nearest' )

   # plot viterbi traceback
   tr_len = len(TR[0])-1

   # draw the viterbi trace
   if tr_len > 2:
      ax1.scatter( TR[0][0], TR[1][0], c='k' )
      ax1.scatter( TR[0][tr_len], TR[1][tr_len], c='k' )

      ax1.plot( TR[0], TR[1], 'k-')

      # Create a Rectangle patch for Viterbi window
      # rect = patches.Rectangle((50,100),40,30,linewidth=1,edgecolor='r',facecolor=None)
      # ax1.add_patch(rect)

      # Add legend
      cbar = ax1.figure.colorbar(im, ax=ax1)

   plt.show()
   return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# Parse args
if (len(sys.argv) == 1):
   print("Usage: <tsv_matrix_1>");
   sys.exit(0);
else:
   for i in range(1, len(sys.argv)):
      arg = sys.argv[i];
      # apply options
      if (arg.startswith("-")):
         if arg in opts.keys():
            opts[arg] = True;
      else:
         tsv_matrix.append(sys.argv[i]);

print(opts)
print(tsv_matrix)

# number of matrices
N = len(tsv_matrix);

# Load matrix (normal and special)
NM_MX = [];
SP_MX = [];
TRACE = [];
for i in range(N):
   N_MX, S_MX, TR = load_matrix(tsv_matrix[i]);
   NM_MX.append(N_MX);
   SP_MX.append(S_MX);
   TRACE.append(TR);

# Find min, max, range of matrix
min_val = np.inf;
max_val = -np.inf;
for i in range(N):
   mx_min, mx_max = minmax_matrix(NM_MX[i], SP_MX[i]);
   min_val = min( min_val, mx_min );
   max_val = max( max_val, mx_max );
rang = max_val - min_val;
print("min:", min_val, "max:", max_val, "range:", rang);

# Remove infinite values from matrix and replace with (min/max values * multiplier)
for i in range(N):
   remove_inf_matrix(NM_MX[i], SP_MX[i], min_val, max_val);
   next;

# if diff, then take the difference of the two matrices
if opts["-diff"]:
   if len(tsv_matrix) == 2:
      NM_MX[0] = diff_matrix( NM_MX[0], NM_MX[1] );
      SP_MX[0] = diff_matrix( SP_MX[0], SP_MX[1] );
      NM_MX.pop(1)
      SP_MX.pop(1)
      N -= 1
   else:
      print("To use -diff, requires two matrices")
      exit(0)

# if add, then take the sum of the two matrices
if opts["-add"]:
   if len(tsv_matrix) == 2:
      NM_MX[0] = add_matrix( NM_MX[0], NM_MX[1] );
      SP_MX[0] = add_matrix( SP_MX[0], SP_MX[1] );
      NM_MX.pop(1)
      SP_MX.pop(1)
      N -= 1
   else:
      print("To use -diff, requires two matrices")
      exit(0)

# if eq, then compare equality of two matrices
if opts["-eq"]:
   if len(tsv_matrix) == 2:
      NM_MX[0] = eq_matrix( NM_MX[0], NM_MX[1] );
      SP_MX[0] = eq_matrix( SP_MX[0], SP_MX[1] );
      NM_MX.pop(1)
      SP_MX.pop(1)
      N -= 1
   else:
      print("To use -eq, requires two matrices")
      exit(0)

# title
title1 = "MATCH\n(min=%f, max=%f)" % (min_val, max_val)
title2 = title1

# # Normalize matrix
# for i in range(N):
#    normalize_matrix(NM_MX[i], min_val, max_val, rang)
#    normalize_matrix(SP_MX[i], min_val, max_val, rang)
#    next

# # Normalize matrix by antidiag (for each channel)
# for i in range(N):
#    print("normalizing antidiags...")
#    normalize_antidiags(NM_MX[i])
#    next

# Transpose matrix so that query is on x-axis, model is on y-axis
for i in range(N):
   # NM_MX[i] = NM_MX[i].T
   next

# Split matrix into state matrices: match, delete, insert
MAT_MX = []
INS_MX = []
DEL_MX = []
for i in range(N):
   MAT_MX.append( NM_MX[i][:,:,0] )
   INS_MX.append( NM_MX[i][:,:,1] )
   DEL_MX.append( NM_MX[i][:,:,2] )
   next

print("shape", MAT_MX[0].shape)

# Print matrices
for i in range(N):
   # print( 'MAT:\n', MAT_MX[i] )
   # print( 'INS:\n', INS_MX[i] )
   # print( 'DEL:\n', DEL_MX[i] )
   # print( 'SPECIAL:\n', SP_MX[i] )
   next

# Output matrices as heatmaps
for i in range(N):
   if opts["-4"]:
      output_heatmap(title1, MAT_MX[i], INS_MX[i], DEL_MX[i], SP_MX[i], min_val, max_val)
   if opts["-tr"]:
      output_heatmap_trace(title2, MAT_MX[i], TRACE[i], min_val, max_val)
   next
