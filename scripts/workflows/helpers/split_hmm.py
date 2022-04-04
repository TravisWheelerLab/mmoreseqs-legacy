#!/usr/bin/env python
###############################################################################
#  - FILE: 	split_hmm.py
#  - DESC:  	Split .hmm into subfiles of given size.
###############################################################################

import os
import sys

# parse commandline args
pwd = os.getcwd()
eval = ""
if len(sys.argv) == 2:
    in_file = sys.argv[i]
else:
    print("Usage: <hmm_file>")
    sys.exit(1)

in_fp = open(in_file, "r")

num_hmm = 0
out_file = "{}.{}".format(in_file, num_hmm)
out_fp = open(out_file, "w")

for line in in_fp:
    out_fp.write(line)
    if line == "//":
        num_hmm += 1
        out_fp.close()
        out_file = "{}.{}".format(in_file, num_hmm)
        out_fp = open(out_file, "w")
