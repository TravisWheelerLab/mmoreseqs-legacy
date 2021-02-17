#!/usr/bin/env python
###############################################################################
#    FILE: 	split_hmm.py
#   BRIEF: 	Split .hmm by batch size
#
#  AUTHOR: 	Dave Rich
###############################################################################

import os,sys

# parse commandline args
pwd = os.getcwd()
eval = ""
if len(sys.argv) == 3:
	in_file = sys.argv[1]
	batch_size = int(sys.argv[2])
else:
   print("Usage: <hmm_file> <batch_size>")
   sys.exit()

in_fp = open(in_file, "r")

# extract file name and extension
file_parts = in_file.split(".")
file_ext = file_parts[len(file_parts)-1]
file_name = ""
file_parts = file_parts[0:len(file_parts)-1]
for part in file_parts:
	file_name += part

print("# FILE: {}".format(in_file))
print("# FILE_NAME: {}".format(file_name))
print("# FILE_EXT: {}".format(file_ext))
print("# BATCH_SIZE: {}".format(batch_size))

num_models = 1
num_models_in_batch = 1
num_batches = 1
out_file = "{}.{}.{}".format(file_name, num_batches-1, file_ext)
out_fp = open(out_file, "w")

for line in in_fp:
	out_fp.write(line)
	if line.startswith("//"):
		# when number of models in current batch reaches size threshold open new file
		if num_models_in_batch >= batch_size:
			out_fp.close()
			num_batches += 1
			out_file = "{}.{}.{}".format(file_name, num_batches-1, file_ext)
			out_fp = open(out_file, "w")	
			num_models_in_batch = 0
		num_models_in_batch += 1
		num_models += 1
out_fp.close()

print("# TOTAL_MODELS: {}".format(num_models))
print("# TOTAL_BATCHES: {}".format(num_batches))
print("# [done].")
