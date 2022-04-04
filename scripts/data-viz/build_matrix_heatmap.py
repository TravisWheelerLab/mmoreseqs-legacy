#!/usr/bin/env python
###############################################################################
#  - FILE: prob_to_bitmap
#  - DESC:  Takes in tsv of dp matrix and creates a bitmap.
###############################################################################

import sys
import argparse
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cmaps
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn import preprocessing

# set resolution of output
mpl.rcParams['figure.dpi'] = 200

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
    "-1": False,         # include only Match matrices
    "-dot": False,       # dot matrix
    "-4": False,         # include all 4 matrices
    "-tr": False,        # include traceback
    "-diff": False,      # diff of two matrices
    "-sum": False,       # sum of two matrices
    "-eq": False,
    "-title": False,     # set title
    "-S": False,         # save figure
    "-ND": False         # no display
}

tsv_matrix = []

tname = "Target Sequence"
qname = "Query Sequence"
minmax = [0, 0]

# If value is infinity, replace it with another value


def replace_inf(val, replace):
    if val == np.inf:
        return replace
    if val == -np.inf:
        return -replace
    return val

# Load matrix in from .tsv file


def load_matrix(tsv_matrix):
    trace = [[], []]

    TR = None
    main_mx = np.zeros((1, 1, 1))
    spec_mx = np.zeros((1, 1))

    # read every line in file
    with open(tsv_matrix) as f:
        line = f.readline()
        while line:

            # ignore comment lines
            if line[0] == "#":
                pass

            # dimensions of matrix (single line)
            elif line.startswith("X_DIM"):
                line = line.split("\t")
                m = int(line[1])
                n = int(line[2])
                main_mx = np.zeros((n+1, m+1, nm_st))
                print("main_mx:", main_mx.shape)
                spec_mx = np.zeros((m+1, sp_st))
                print("spec_mx:", spec_mx.shape)
                pass

            # traceback indices
            elif line.startswith("X_trace"):
                line = f.readline()
                while line and not line.startswith("/"):
                    line = line.split("\t")
                    trace[0].append(int(line[1]))
                    trace[1].append(int(line[2]))

                    line = f.readline()
                pass

            # matrix values
            elif line.startswith("X_MATRIX"):
                line = f.readline()
                lines_read = 0

                # read into matrix until EOF or X keyword
                while line and not line.startswith("X"):
                    # label which matrix it goes to
                    fields = line.split()

                    # ignore if line has improper number of fields
                    if line.startswith("#") or (len(fields) < n + 2):
                        # print("# Ignoring line...")
                        # print("# NUM_FIELDS:", len(fields))
                        line = f.readline()
                        continue

                    # label tells whether {M,I,D}
                    label = fields[0]

                    # read in row from main matrix
                    if label in nm_dict.keys():
                        lines_read += 1

                        mx_cur = main_mx
                        st_cur = nm_dict[label]
                        row_cur = int(fields[1])
                        # for cell in row
                        for i in range(n+1):
                            try:
                                mx_cur[i, row_cur, st_cur] = (
                                    float)(fields[i+2])
                            except Exception as e:
                                print("ERROR:", e)
                                print("mx_cur", i, row_cur, st_cur)
                                print("main_mx:", mx_cur.shape)
                                exit(1)
                        pass

                    line = f.readline()
                pass

            elif line.startswith("X_SPECIAL"):
                line = f.readline()
                while line and not line.startswith("/"):
                    # label which row
                    line = line.split()
                    label = line[0]

                    # read in row from special state matrix
                    if label in sp_dict.keys():
                        mx_cur = spec_mx
                        st_cur = sp_dict[label]

                        for i in range(m+1):
                            mx_cur[i][st_cur] = (float)(line[i+1])
                        pass

                    line = f.readline()
                pass

            line = f.readline()

    return main_mx, spec_mx.T, trace

# Take cell-wise difference of two matrices


def diff_matrix(mat1, mat2):
    res = np.zeros_like(mat1)

    for i, x in np.ndenumerate(mat1):
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

    for i, x in np.ndenumerate(mat1):
        if (mat1[i] == -mat2[i]):
            res[i] = 0
        else:
            res[i] = mat1[i] + mat2[i]

    return res

# Test equality of two matrices (within tolerance)


def eq_matrix(mat1, mat2, tol=0.01):
    res = np.zeros_like(mat1)

    for i, x in np.ndenumerate(mat1):
        if abs(mat1[i] - mat2[i]) < tol:
            res[i] = 0
        else:
            res[i] = 1

    return res

# Find non-infinite min/max


def minmax_matrix(main_mx, spec_mx):
    min_val = np.inf
    max_val = -np.inf

    # find min and max non-infinite values
    for i, x in np.ndenumerate(main_mx):
        val = main_mx[i]
        if (val != np.inf and val != -np.inf):
            if min_val > val:
                min_val = val
            if max_val < val:
                max_val = val

    # find min and max non-infinite values
    for i, x in np.ndenumerate(spec_mx):
        val = spec_mx[i]
        if (val != np.inf and val != -np.inf):
            if min_val > val:
                min_val = val
            if max_val < val:
                max_val = val

    return min_val, max_val

# Replace infinite values with non-infinite values: inf=max*mult, -inf=min*mult


def remove_inf_matrix(main_mx, spec_mx, min_val, max_val):
    multiplier = 1.2
    # remove infinite values
    for i, x in np.ndenumerate(main_mx):
        val = main_mx[i]
        if (val == np.inf):
            main_mx[i] = max_val * multiplier
        elif (val == -np.inf):
            main_mx[i] = min_val * multiplier

    # remove infinite values
    for i, x in np.ndenumerate(spec_mx):
        val = spec_mx[i]
        if (val == np.inf):
            spec_mx[i] = max_val * multiplier
        elif (val == -np.inf):
            spec_mx[i] = min_val * multiplier
    return None

# Normalize matrix (given min, max, and range)


def normalize_matrix(main_mx, min_val, max_val, rang):
    main_mx -= min_val
    main_mx *= (1 / rang)
    main_mx = np.exp(main_mx)
    # main_mx *= 255
    # main_mx = main_mx.astype(int)
    return None

# Normalize all anti-diagonals of matrix


def normalize_antidiags(main_mx):
    print(main_mx.shape)
    x, y, C = main_mx.shape
    min_corner = min(x, y)
    max_corner = max(x, y)
    num_diags = x+y-1

    # for each state (M,I,D)
    for c in range(C):
        num_cells = 0
        start_i = 0
        M = main_mx

        print("num_diags:", num_diags)

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
                    if M[i, j, c] == np.inf or M[i, j, c] == -np.inf:
                        next
                    if M[i, j, c] > max_val:
                        max_val = M[i, j, c]
                    if M[i, j, c] < min_val:
                        min_val = M[i, j, c]

            # replace inf vals with min/max vals
            for i in range(start_i, start_i+num_cells):
                for c in range(C):
                    j = d-i
                    # M[i,j,c] = d
                    cell_cnt += 1
                    # print(d,i,j)
                    if M[i, j, c] == np.inf:
                        M[i, j, c] = max_val
                    if M[i, j, c] == -np.inf:
                        M[i, j, c] = min_val

            range_val = max_val - min_val
            # print("PRE d:",d,"range_val:",range_val,"min_val:",min_val, "max_val:", max_val)

            # normalize the diagonal
            for i in range(start_i, start_i+num_cells):
                for c in range(C):
                    j = d-i
                    if (range_val != 0):
                        M[i, j, c] = ((M[i, j, c] - min_val) / range_val)
                    else:
                        M[i, j, c] = 0.0
            # print("POST d:",d,"range_val:",range_val,"min_val:",min_val, "max_val:", max_val)
    return None

# Output matrix at heatmap (1 grid, MATCH)


def output_dotmap_1(title, mat_mx, ins_mx, del_mx, spec_mx, vmin, vmax, file="", default="Match State Bitscore Heatmap"):
    # matplotlib defaults
    # mpl.rcParams['figure.dpi'] = 300

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 6}

    mpl.rc('font', **font)

    # get colormap
    my_cmap_name = 'jet'
    my_cmap = cmaps.get_cmap(my_cmap_name)

    # dot sizes
    dot_range = [1, 20]
    # resolution of
    res = 3
    whitepad = 5
    my_dim = mat_mx.shape

    # layout plot
    color_mx = mat_mx.copy()
    c_minmax = [color_mx.min(), color_mx.max()]
    c_range = c_minmax[1] - c_minmax[0]

    fig, ((ax1)) = plt.subplots(1, 1)
    ax1.set_title(title)
    ax1.set_ylabel(qname)
    ax1.set_xlabel(tname)

    # log function
    # color_mx = color_mx[:] - color_mx.min()
    # color_mx = np.exp2(color_mx)

    # # normalized nonlinear step function
    # bps = [c_minmax[0], 20, 50, c_minmax[1]]   # breakpoints in data
    # prcs = [0, 0.05, 0.50, 1.00]            # percentages of colormap usage
    # color_mx[ color_mx == bps[0] ] = 0.0
    # # subtot_prc = 0.0
    # print(color_mx)
    # for i in range(1, len(bps)):
    #    inner_prc  = prcs[i] - prcs[i-1]
    #    inner_rng = bps[i] - bps[i-1]
    #    # inner_mx = color_mx[ (color_mx > bps[i-1]) & (color_mx <= bps[i])  ]
    #    # inner_mx = (((inner_mx - bps[i-1]) / inner_rng) * inner_prc) + prcs[i]
    #    color_mx[ (color_mx > bps[i-1]) & (color_mx <= bps[i])  ] = (((color_mx[ (color_mx > bps[i-1]) & (color_mx <= bps[i])  ] - bps[i-1]) / inner_rng) * inner_prc) + prcs[i-1]
    # print(color_mx)

    # color_mx[ color_mx > bps[i-1] and color_mx <= bps[i]  ] = ((( color_mx[ color_mx <= 20 ] - c_minmax[0] ) / (20 - c_minmax[0] )) * 0.05) + 0.0
    # color_mx[ color_mx > 20 ]  = ((( color_mx[ color_mx > 20 ] - 20 ) / (c_minmax[1] - 20)) * 0.95) + 0.05

    # output heatmap
    img = ax1.imshow(color_mx, cmap=my_cmap_name, interpolation=None,
                     alpha=((color_mx - c_minmax[0])/c_range))

    # colored level contours
    my_y = range(my_dim[0])
    my_x = range(my_dim[1])
    my_range = minmax[1] - minmax[0]
    # my_levels = np.linspace(minmax[0], minmax[1], 21)
    # ax1.contour( my_x, my_y, mat_mx, cmap=my_cmap_name, levels=20, linestyles='-')

    # # find dot minmax
    # dot_minmax = [np.inf, -np.inf]
    # for i in my_x[::res]:
    #    my_scores = mat_mx[:,i][::res]
    #    dot_minmax[0] = min( dot_minmax[0], my_scores.min() )
    #    dot_minmax[1] = max( dot_minmax[1], my_scores.max() )
    # print("dot_minmax:", dot_minmax[0], dot_minmax[1])
    # dot_range = dot_minmax[1] - dot_minmax[0]

    # # create spaced dots
    # dotmax = -np.inf
    # for i in my_x[::res]:
    #    my_scores = mat_mx[:,i][::res]
    #    my_colors = my_scores[:]
    #    my_sizes  = (((my_colors[:] - dot_minmax[0])/dot_range) * 2.5)
    #    ax1.scatter( [i] * len(my_y[::res]), my_y[::res], cmap=my_cmap_name, c=mat_mx[:,i][::res], marker='o', s=my_sizes, vmin=dot_minmax[0], vmax=dot_minmax[1])

    # add ticks
    # dim = max(mat_mx.shape[0])
    # yticks = math.floor(math.log())
    ax1.xaxis.grid(True, which='minor')
    plt.grid(color='gray', linestyle='-', linewidth=1, alpha=0.25)

    # Add legend
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(img, cax=cax)
    # cbar = ax1.figure.colorbar(img, ax=ax1)

    plt.tight_layout()
    if (opts["-S"]):
        dest = "{}".format(name)
        dest = "heatmap.jpg"
        plt.savefig(dest)
        print("Figure saved to '{}'...".format(dest))

    if not (opts["-ND"]):
        print("Displaying figure...")
        plt.show()

    return None

# Output matrix at heatmap (1 grid, MATCH)


def output_heatmap_1(title, mat_mx, ins_mx, del_mx, spec_mx, vmin, vmax, file="", default="Match State Bitscore Heatmap"):
    # matplotlib defaults
    # mpl.rcParams['figure.dpi'] = 300

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 6}

    mpl.rc('font', **font)

    # layout plot
    fig, ((ax1)) = plt.subplots(1, 1)
    ax1.set_title(title)
    ax1.set_ylabel(qname)
    ax1.set_xlabel(tname)
    img = ax1.imshow(mat_mx, cmap='jet', interpolation='nearest',
                     vmin=minmax[0], vmax=minmax[1])

    # add ticks
    # dim = max(mat_mx.shape[0])
    # yticks = math.floor(math.log())
    ax1.xaxis.grid(True, which='minor')
    plt.grid(color='gray', linestyle='-', linewidth=1, alpha=0.25)

    # Add legend
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(img, cax=cax)
    # cbar = ax1.figure.colorbar(img, ax=ax1)

    plt.tight_layout()
    if (opts["-S"]):
        dest = "{}".format(name)
        dest = "heatmap.jpg"
        plt.savefig(dest)
        print("Figure saved to '{}'...".format(dest))

    if not (opts["-ND"]):
        print("Displaying figure...")
        plt.show()

    return None

# Output matrix at heatmap


def output_heatmap_4(title, mat_mx, ins_mx, del_mx, spec_mx, vmin, vmax, file="", default="test"):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.set_title(title)
    ax1.imshow(mat_mx, cmap='jet', interpolation='nearest')
    ax2.set_title("INSERT")
    ax2.imshow(ins_mx, cmap='jet', interpolation='nearest')
    ax3.set_title("DELETE")
    ax3.imshow(del_mx, cmap='jet', interpolation='nearest')
    ax4.set_title("SPECIAL")
    im = ax4.imshow(spec_mx, cmap='jet', interpolation='nearest')

    # Add legend
    cbar = ax1.figure.colorbar(im, ax=ax4)

    plt.tight_layout()
    if (opts["-S"]):
        dest = "{}".format(name)
        dest = "heatmap.jpg"
        plt.savefig(dest)
        print("Figure saved to '{}'...".format(dest))

    if not (opts["-ND"]):
        print("Displaying figure...")
        plt.show()

    return None

# Output heatmap of single matrix with traceback on top


def output_heatmap_trace(title, mat_mx, TR, vmin, vmax, file="", name="test"):
    fig, ax1 = plt.subplots(1, 1)

    # plot heatmap
    ax1.set_title(title)
    im = ax1.imshow(mat_mx, cmap='jet', interpolation='nearest')

    # plot viterbi traceback
    if (TR != None):
        tr_len = len(TR[0])-1

    # draw the viterbi trace
    if TR != None and tr_len >= 2:
        print("making trace...\n")
        ax1.scatter(TR[0][0], TR[1][0], c='white', s=3)
        ax1.scatter(TR[0][tr_len], TR[1][tr_len], c='white', s=3)

        #ax1.plot( TR[0], TR[1], linestyle='-', linewidth=1, color='black')
        ax1.plot(TR[0], TR[1], linestyle='-', linewidth=2, color='white')
        ax1.plot(TR[0], TR[1], linestyle='-', linewidth=1, color='black')

        # Create a Rectangle patch for Viterbi window
        # rect = patches.Rectangle((50,100),40,30,linewidth=1,edgecolor='r',facecolor=None)
        # ax1.add_patch(rect)

        # Add legend
        #cbar = ax1.figure.colorbar(im, ax=ax1)

    if (opts["-S"]):
        dest = "{}".format(name)
        dest = "heatmap.jpg"
        plt.savefig(dest)
        print("Figure saved to '{}'...".format(dest))

    if not (opts["-ND"]):
        print("Displaying figure...")
        plt.show()

    return None


##############################################################################
###########################         MAIN         #############################
##############################################################################

# default location to save files
# file = "/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/data-vis/heatmaps/"
file = ""
title = "heatmap.viz"

# Parse args
if (len(sys.argv) == 1):
    print("Usage: <tsv_matrix_1>")
    sys.exit(0)
else:
    i = 1
    while (i < len(sys.argv)):
        arg = sys.argv[i]
        # apply options
        if (arg == "-title"):
            i += 1
            title = sys.argv[i]
        if (arg == "-tname"):
            i += 1
            tname = sys.argv[i]
        if (arg == "-qname"):
            i += 1
            qname = sys.argv[i]
        if (arg == "-range"):
            i += 1
            minmax[0] = float(sys.argv[i])
            i += 1
            minmax[1] = float(sys.argv[i])
        elif (arg.startswith("-")):
            if arg in opts.keys():
                opts[arg] = True
        else:
            tsv_matrix.append(sys.argv[i])
        i += 1

print("OPTS:", opts)
print("MATRICES:", tsv_matrix)

# number of matrices
N = len(tsv_matrix)

# Load matrix (normal and special)
main_mx = []
spec_mx = []
trace = []
for i in range(N):
    N_MX, S_MX, TR = load_matrix(tsv_matrix[i])
    main_mx.append(N_MX)
    spec_mx.append(S_MX)
    trace.append(TR)

# Find min, max, range of matrix
min_val = np.inf
max_val = -np.inf
for i in range(N):
    mx_min, mx_max = minmax_matrix(main_mx[i], spec_mx[i])
    min_val = min(min_val, mx_min)
    max_val = max(max_val, mx_max)

rang = max_val - min_val
print("min:", min_val, "max:", max_val, "range:", rang)
minmax[0] = min(minmax[0], min_val)
minmax[1] = max(minmax[1], max_val)

# Remove infinite values from matrix and replace with (min/max values * multiplier)
for i in range(N):
    remove_inf_matrix(main_mx[i], spec_mx[i], min_val, max_val)
    next

# if diff, then take the difference of the two matrices
if opts["-diff"]:
    if len(tsv_matrix) == 2:
        main_mx[0] = diff_matrix(main_mx[0], main_mx[1])
        spec_mx[0] = diff_matrix(spec_mx[0], spec_mx[1])
        main_mx.pop(1)
        spec_mx.pop(1)
        N -= 1
    else:
        print("To use -diff, requires two matrices")
        exit()

# if add, then take the sum of the two matrices
if opts["-sum"]:
    if len(tsv_matrix) == 2:
        main_mx[0] = add_matrix(main_mx[0], main_mx[1])
        spec_mx[0] = add_matrix(spec_mx[0], spec_mx[1])
        main_mx.pop(1)
        spec_mx.pop(1)
        N -= 1
    else:
        print("To use -diff, requires two matrices")
        exit()

# if eq, then compare equality of two matrices
if opts["-eq"]:
    if len(tsv_matrix) == 2:
        main_mx[0] = eq_matrix(main_mx[0], main_mx[1])
        spec_mx[0] = eq_matrix(spec_mx[0], spec_mx[1])
        main_mx.pop(1)
        spec_mx.pop(1)
        N -= 1
    else:
        print("To use -eq, requires two matrices")
        exit()

# # Normalize matrix
# for i in range(N):
#    normalize_matrix(main_mx[i], min_val, max_val, rang)
#    normalize_matrix(spec_mx[i], min_val, max_val, rang)
#    next

# # Normalize matrix by antidiag (for each channel)
# for i in range(N):
#    print("normalizing antidiags...")
#    normalize_antidiags(main_mx[i])
#    next

# Transpose matrix so that sequence(target) is on x-axis, model(query) is on y-axis
for i in range(N):
    # main_mx[i] = main_mx[i].T
    next

# Split matrix into state matrices: match, delete, insert
mat_mx = []
ins_mx = []
del_mx = []
for i in range(N):
    mat_mx.append(main_mx[i][:, :, 0])
    ins_mx.append(main_mx[i][:, :, 1])
    del_mx.append(main_mx[i][:, :, 2])
    next

print("shape", mat_mx[0].shape)

# Print matrices
for i in range(N):
    # print( 'MAT:\n', mat_mx[i] )
    # print( 'INS:\n', ins_mx[i] )
    # print( 'DEL:\n', del_mx[i] )
    # print( 'SPECIAL:\n', spec_mx[i] )
    next

# Output matrices as heatmaps
for i in range(N):
    name = tsv_matrix[i].split("/")
    name = name[len(name)-1]

    if opts["-1"]:
        output_heatmap_1(
            title, mat_mx[i], ins_mx[i], del_mx[i], spec_mx[i], min_val, max_val, file, name)
    if opts["-dot"]:
        output_dotmap_1(
            title, mat_mx[i], ins_mx[i], del_mx[i], spec_mx[i], min_val, max_val, file, name)
    if opts["-4"]:
        output_heatmap_4(
            title, mat_mx[i], ins_mx[i], del_mx[i], spec_mx[i], min_val, max_val, file, name)
    if opts["-tr"]:
        output_heatmap_trace(
            title, mat_mx[i], trace[i], min_val, max_val, file, name)
    next
