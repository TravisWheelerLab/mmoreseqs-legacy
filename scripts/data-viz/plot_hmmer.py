#!/usr/bin/env python
###############################################################################
#  - FILE: filter_results.py
#  - DESC:  Filter results
###############################################################################

from __future__ import print_function
import sys
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# write to stderr


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    return

# load hmmer


def load_hmmer_time(filename):
    # line counter
    l_cnt = 0
    # field counter
    f_cnt = 0
    # output data
    data = []

    with open(filename, "r") as fp:
        for line in fp:
            # check for comment lines or empty lines
            if (line.startswith("##")):
                continue

            fields = line.split()
            # check that there are a proper number of fields
            if (len(fields) != 7):
                continue

            total_cells = (int(fields[0])+1) * (int(fields[2])+1)
            fields.append(total_cells)

            data.append(fields)
        return data

# append result_id, target_id, query_id and write out


def append_hmmer_time(filename, res_lookup, t_lookup, q_lookup):
    # line counter
    l_cnt = 0
    # field counter
    f_cnt = 0

    with open(filename, "r") as fp:
        for line in fp:
            # check for comment lines or empty lines
            if (line.startswith("##")):
                continue
            # break up line into space/tab separated components
            line = line.rstrip()
            fields = line.split()
            # check that there are a proper number of fields
            if (len(fields) != 7):
                continue

            try:
                t_name = fields[1]
                q_name = fields[3]

                t_id = t_lookup[t_name]
                q_id = q_lookup[q_name]

                # only report hmmer result if it is also in mmseqs results
                if (t_id, q_id) in res_lookup.keys():
                    res_id = res_lookup[(t_id, q_id)]

                    print(f'{res_id} {t_id} {q_id} {line}')

                l_cnt += 1
                if (l_cnt % 100000 == 0):
                    eprint(f'# {l_cnt} entries processed...')

            except:
                print(f'# ERROR occurred on line {l_cnt}: "{line}" ')
                eprint(f'# ERROR occurred on line {l_cnt}: "{line}" ')
            continue
    return


# load lookup dict
def load_lookup(filename):
    # line counter
    l_cnt = 0
    # field counter
    f_cnt = 0
    # output lookup dict
    name_lookup = {}
    id_lookup = {}

    # read file line-by-line
    with open(filename, "r") as fp:
        for line in fp:
            # check for comment or improper format lines
            if (line.startswith("#")):
                continue
            fields = line.split()
            if (len(fields) != 3):
                continue

            id_ = fields[0]
            offset = fields[1]
            name = fields[2]
            name_lookup[id_] = name
            id_lookup[name] = id_
            pass
        pass
    pass
    return name_lookup, id_lookup

# m8+ results lookup dict


def load_result_lookup(filename):
    # line counter
    l_cnt = 0
    # field counter
    f_cnt = 0
    # output lookup dict
    res_lookup = {}

    # read file line-by-line
    with open(filename, "r") as fp:
        for line in fp:
            # check for comment or improper format lines
            if (line.startswith("#")):
                continue
            fields = line.split()
            if (len(fields) != 17):
                continue

            res_id = fields[0]
            t_id = fields[1]
            q_id = fields[2]
            res_lookup[(t_id, q_id)] = res_id
        pass
    pass
    return res_lookup


##############################################################################
###########################         MAIN         #############################
##############################################################################

# get args
if ((len(sys.argv) == 2) and (sys.argv[1] == "--defaults")):
    hmmer_time_fname = "/home/devreckas/Data/hmmer-test-stdout.txt"
    target_lookup_fname = "/home/devreckas/Data/profmark-benchmark/temp/cloud/pmark.hmm.idx"
    query_lookup_fname = "/home/devreckas/Data/profmark-benchmark/temp/cloud/pmark.fa.idx"
    results_m8_fname = "/home/devreckas/Data/profmark-benchmark/temp/cloud/mmseqs_res.m8+"
elif (len(sys.argv) != 5):
    print("Usage: <hmmer_times> <target_lookup> <query_lookup> <results_m8+>")
    exit(0)
else:
    hmmer_time_fname = sys.argv[1]
    target_lookup_fname = sys.argv[2]
    query_lookup_fname = sys.argv[3]
    results_m8_fname = sys.argv[4]
pass

# test print
# print("This println goes to stdout.")
# eprint("This print goes to stderr.")
# exit(0)

# load data
target_name_lookup, target_id_lookup = load_lookup(target_lookup_fname)
query_name_lookup, query_id_lookup = load_lookup(query_lookup_fname)
res_id_lookup = load_result_lookup(results_m8_fname)
# hmmer_data = load_hmmer_time(hmmer_time_fname)

append_hmmer_time(hmmer_time_fname, res_id_lookup,
                  target_id_lookup, query_id_lookup)

# # linear spacing
# x = range( len(mmseqs_plus_data) )

# # sort data by size
# hmmer_data.sort( key=lambda r:int(r[7]) )

# # plot time data against problem size
# fwdbck_time = []
# size = []
# for data in hmmer_data:
# 	fwdbck_time.append( float(data[4]) )
# 	size.append( int(data[7]) )

# plt.scatter( size, fwdbck_time, marker='o', s=2 )
# plt.show()
