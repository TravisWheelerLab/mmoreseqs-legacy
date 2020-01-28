 #!/usr/bin/env python

 ###############################################################################
 #  @file  benchmark_linegraph
 #  @brief Render benchmark scores as a line graph
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

filename = ""
if len(sys.argv) == 2:
   filename = sys.argv[1]
else:
   print("Usage: <results_filename>")
   sys.exit(0)

name       = ""
viterbi    = []
fwd        = []
bck        = []
cloud_fwd  = []
cloud_bck  = []
alpha      = []
beta       = []
perc_cells = []
perc_wind  = []
perc_tot   = []

# read in file
line_cnt = 0
with open(filename, "r") as fp:
   for line in fp:
      line = line.split()
      print(line)
      name = line[0]
      viterbi.append( float(line[1]) )
      fwd.append( float(line[2]) )
      bck.append( float(line[3]) )
      cloud_fwd.append( float(line[4]) )
      cloud_bck.append( float(line[5]) )
      alpha.append( float(line[6]) )
      beta.append( int(line[7]) )
      perc_cells.append( float(line[8]) * 100 )
      perc_wind.append( float(line[9]) * 100 )
      perc_tot.append( 100 )

      line_cnt += 1
      print(line_cnt)


# render results
color_set = []
line_styles = []
plt.plot(alpha, viterbi, 'r--', label="viterbi score")
# plt.plot(alpha, alpha, 'b--', label="alpha score")
plt.plot(alpha, fwd, 'b.', label="forward score")
plt.plot(alpha, fwd, 'b--', label="backward score")
plt.plot(alpha, cloud_fwd, 'g.', label="cloud-fwd score")
plt.plot(alpha, cloud_bck, 'g--', label="cloud-bck score")
plt.plot(alpha, perc_cells, 'k--', label="percent cells computed")
plt.plot(alpha, perc_tot, 'k-', label="total percent")
plt.xlabel('alpha pruning threshold')
plt.show()


