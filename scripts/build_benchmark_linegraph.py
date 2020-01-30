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

filenames = []
eval = ""
if len(sys.argv) >= 2:
   filenames.append(sys.argv[1])
else:
   print("Usage: <results_filename1> <results_filename2> ...")
   sys.exit(0)

if len(sys.argv) >= 3:
   e_val = int(sys.argv[2])
   qlen = int(sys.argv[3])
   tlen = int(sys.argv[4])

# read in file
for filename in filenames:
   hmm_name   = []
   fasta_name = []
   fa_name    = []
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

   line_cnt = 0
   with open(filename, "r") as fp:
      for line in fp:
         line = line.split()
         print("line:", line)
         name = line[0].split("/")
         name = name[len(name)-1]
         hmm_name.append(name)
         name = line[1].split("/")
         name = name[len(name)-1]
         fasta_name.append(name)
         viterbi.append( float(line[2]) )
         fwd.append( float(line[3]) )
         bck.append( float(line[4]) )
         cloud_fwd.append( float(line[5]) )
         cloud_bck.append( float(line[6]) )
         alpha.append( float(line[7]) )
         beta.append( int(line[8]) )
         perc_cells.append( float(line[9]) * 100 )
         perc_wind.append( float(line[10]) * 100 )
         perc_tot.append( 100 )

         line_cnt += 1
         print(line_cnt)

   # render results
   plt.subplot(2,1,1)
   title = "{} || {}".format(hmm_name[0], fasta_name[0])
   plt.title(title)
   plt.plot(alpha, viterbi, 'r--', label="viterbi score")
   # plt.plot(alpha, alpha, 'b--', label="alpha score")
   plt.plot(alpha, fwd, 'b.', label="forward score")
   plt.plot(alpha, fwd, 'b--', label="backward score")
   plt.plot(alpha, cloud_fwd, 'g.', label="cloud-fwd score")
   plt.plot(alpha, cloud_bck, 'g--', label="cloud-bck score")
   plt.ylabel("score")

   plt.subplot(2,1,2)
   plt.plot(alpha, perc_cells, 'k--', label="percent cells computed")
   plt.plot(alpha, perc_tot, 'k-', label="total percent")
   plt.ylabel("percent of cells computed")

   plt.xlabel('alpha pruning threshold')
   plt.savefig("data-vis/{}.jpg".format(fasta_name[0]))
   # plt.show()


