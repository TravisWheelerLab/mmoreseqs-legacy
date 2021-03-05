#!/usr/bin/env python
###############################################################################
#    FILE: 	split_fasta.py
#   BRIEF: 	Split .fasta file
#
#  AUTHOR: 	Dave Rich
###############################################################################

import os,sys

# parse commandline args
pwd = os.getcwd()
eval = ""

if len(sys.argv) < 3:
   print("Usage: <fasta_file> <batch_size>")
   sys.exit(1)
if len(sys.argv) == 3:
	in_file = sys.argv[1]
	batch_size = int(sys.argv[2])
if len(sys.argv) == 4:
	out_dir = sys.argv[3]

# get file name and extension
file_parts = in_file.split(".")
file_ext = file_parts[len(file_parts)-1]
file_parts = file_parts[0:len(file_parts)-1]
file_name = ""
for part in file_parts:
	file_name += part + "."

# print("# FILE: {}".format(in_file))
# print("# FILE_NAME: {}".format(file_name))
# print("# FILE_EXT: {}".format(file_ext))
# print("# BATCH_SIZE: {}".format(batch_size))

in_fp = open(in_file, "r")
num_batches = 1
num_models_in_batch = 0
num_models = 0
out_file = "{}{}.{}".format(file_name, file_ext, num_batches-1)
out_fp = open(out_file, "w")

for line in in_fp:
	# if at the beginning of a new fasta sequece
	if line.startswith(">"):
		# if batch limit is reached
		if num_models_in_batch >= batch_size:	
			out_fp.close()
			num_batches += 1
			out_file = "{}{}.{}".format(file_name, file_ext, num_batches-1)
			out_fp = open(out_file, "w")
			num_models_in_batch = 0
		num_models += 1
		num_models_in_batch += 1	

	out_fp.write(line)
in_fp.close()
out_fp.close()

print("# TOTAL_BATCHES: {}".format(num_batches))
print("# TOTAL_MODELS: {}".format(num_models))
print("# [done].")

