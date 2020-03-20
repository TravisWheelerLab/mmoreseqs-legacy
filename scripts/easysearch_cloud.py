 #!/usr/bin/env python

###############################################################################
#  @file prob_to_bitmap
#  @brief Takes in tsv of dp matrix and creates a bitmap.
#
#  @author Dave Rich
#  @bug Lots.
###############################################################################

import os, sys
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# get .fasta in file
def copy_fasta(infile, query, outfile):
   print("copy fasta")
   is_in_fasta = False
   f_out = open(outfile, "w")

   with open(infile, "r") as f_in:
      line = f_in.readline()
      i = 0
      while line:
         if line.startswith(">"):
            if line.find(query) != -1:
               print("QUERY FOUND")
               is_in_fasta = True
               f_out.write(line)
               print(i, line)
               line = f_in.readline()
               i += 1
               continue

         if is_in_fasta == True:
            if line.startswith(">"):
               f_out.close()
               return
            print(i, line)
            f_out.write(line)

         line = f_in.readline()
         i += 1



##############################################################################
###########################         MAIN         #############################
##############################################################################

opts = {}
# Parse args
if (len(sys.argv) == 1):
   print("Usage: <tsv_matrix_1>");
   sys.exit(0)
else:
   i = 1
   while (i < len(sys.argv) ):
      arg = sys.argv[i]
      i += 1
      # apply options
      if (arg.startswith("-")):
         optarg = sys.argv[i]
         opts[arg] = optarg
         i += 1

for key in opts.keys():
   print(key, ":", opts[key])

tmp_folder = "my_tmp"
if not os.path.exists(tmp_folder):
   os.makedirs(tmp_folder)
out_folder = "my_output"
if not os.path.exists(out_folder):
   os.makedirs(out_folder)

results_file = opts["-f"]
query_file   = opts["-q"]
target_file  = opts["-t"]

if not opts["-f"]:
   print("ERROR: No results filename.")

hits = []

print("LOADING all hits...")
# load all hits
with open(results_file, "r") as f:
   line = f.readline()
   line_num = 0
   while line: 

      hit = {}
      line = line.split()

      hit["query"]  = line[0]
      hit["target"] = line[1]

      hit["qstart"] = line[6]
      hit["qend"]   = line[7]

      hit["tstart"] = line[8]
      hit["tend"]   = line[9]

      hits.append(hit)

      line = f.readline()
      line_num += 1

print("len(hits):", len(hits))

print("CREATING fasta/hmm files...")
# create .fasta and .hmm files, if necessary
for i in range(len(hits)):
   hit = hits[i]
   query_name  = hit["query"]
   target_name = hit["target"]

   print(hit["query"], hit["target"])

   # find .fasta query starting line number in file
   # cmd = f"grep -n {hits["query"]} {query_file} | head -1"
   # stdout = os.system(cmd)
   # hits["line_num"] = stdout.split(":")[0]
   # print(hits["query"], hits["line_num"])

   outfile_fa = f"{tmp_folder}/{query_name}.fa"
   if not os.path.exists(outfile_fa):
      print("not found query")
      copy_fasta(query_file, query_name, outfile_fa)

   # create .hmm if doesnt exist
   outfile_fa  = f"{tmp_folder}/{target_name}.fa"
   outfile_hmm = f"{tmp_folder}/{target_name}.hmm"
   if not os.path.exists(outfile_hmm):
      print("not found target")
      copy_fasta(target_file, target_name, outfile_fa)
      cmd = f"hmmbuild {outfile_hmm} {outfile_fa}"
      os.system(cmd)


results_file = f"{out_folder}/my_results.txt"
cloud_search = "cloud_fwdback.exe"
# run cloud search
for i in range(len(hits)):
   hit = hits[i]
   query_name  = hit["query"]
   target_name = hit["target"]

   qstart      = hit["qstart"]
   qend        = hit["qend"]  

   tstart      = hit["tstart"]
   tend        = hit["tend"]  

   cmd = f"./{cloud_search} -P  {target_name}.hmm {query_name}.fa -w {tstart} {qstart} {tend} {qend} -o {results_file} "
   os.system(cmd)



