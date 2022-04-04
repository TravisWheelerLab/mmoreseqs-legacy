#!/bin/usr/python
#####################################################################
# - FILE:  make_dummydb_align.py
# - DESC:   Creates simple dummy MMseqs align database.
# 			 Used in converting result ids from profile results to query results.
#####################################################################

# import os
import sys
# import numpy as np

# identifiers for different database types
my_dbtypes = {}
my_dbtypes["amino"] = ("%c%c%c%c" % (0, 0, 0, 0))
my_dbtypes["nucleotide"] = ("%c%c%c%c" % (1, 0, 0, 0))
my_dbtypes["profile"] = ("%c%c%c%c" % (2, 0, 0, 0))
my_dbtypes["prefilter"] = ("%c%c%c%c" % (7, 0, 0, 0))
my_dbtypes["align"] = ("%c%c%c%c" % (5, 0, 0, 0))
my_dbtypes["generic"] = ("%c%c%c%c" % (12, 0, 0, 0))

#####################################################################
##  MAIN  ###########################################################
#####################################################################

if len(sys.argv) < 3:
    print("./ <db_name> <i:db_size>")
    exit(1)

my_args = {}
my_args["db_type"] = "align"

my_args["db_name"] = sys.argv[1]
my_args["db_size"] = int(sys.argv[2])
if len(sys.argv) > 3:
    my_args["db_type"] = sys.argv[3]

my_files = {}
my_files["db_index"] = "{}.index".format(my_args["db_name"])
my_files["db_0"] = "{}".format(my_args["db_name"])
my_files["db_type"] = "{}.dbtype".format(my_args["db_name"])

fp_index = open(my_files["db_index"], "w")
fp_0 = open(my_files["db_0"], "w")
fp_type = open(my_files["db_type"], "w")

total_offset = 0
entry_size = 0

for i in range(my_args["db_size"]):
    # make entry for database
    # targetID  alnScore  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen
    # [queryOrfStart] [queryOrfEnd] [dbOrfStart] [dbOrfEnd] [alnCigar]
    line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\x00".format(
        i, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    entry_size = len(line)
    fp_0.write(line)

    # make entry for index
    index = "{}\t{}\t{}\n".format(i, total_offset, total_offset)
    total_offset += entry_size
    fp_index.write(index)

# make dbtype file
fp_type.write(my_dbtypes[my_args["db_type"]])

fp_0.close()
fp_index.close()
fp_type.close()

print("# completed successfully.")
