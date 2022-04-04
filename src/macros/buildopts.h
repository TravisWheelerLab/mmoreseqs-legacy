/*******************************************************************************
 *  - FILE:      buildopts.c
 *  - DESC:    Macros which determine build options.
 *******************************************************************************/

/* === BUILD VERSIONING === */
/* Program Name */
#ifndef BUILD_PROGRAM
#define BUILD_PROGRAM "MMORESEQS"
#endif
/* Version Number */
#ifndef BUILD_VERSION
#define BUILD_VERSION "Build Version Unknown"
#endif
/* Name of Version */
#ifndef BUILD_TYPE
#define BUILD_TYPE "Build Type Unknown"
#endif
/* Date of Build */
#ifndef BUILD_DATE
#define BUILD_DATE "March 2021"
#endif
/* Git Hash of Build */
#ifndef BUILD_HASH
#define BUILD_HASH "Build Hash Unknown"
#endif
/* Build Description */
#ifndef BUILD_DESC
#define BUILD_DESC \
  "Heuristic Pruned Forward-Backward for Faster Profile-to-Sequence Search"
#endif
/* Build Repository */
#ifndef BUILD_REPO
#define BUILD_REPO "http://github.com/TravisWheelerLab/fb-pruner/"
#endif
/* Build Copyright */
#ifndef BUILD_COPYRIGHT
#define BUILD_COPYRIGHT \
  "Copyright (C) 2020 Travis Wheeler Lab, University of Montana"
#endif

/* === RESOURCE FOLDER LOCATIONS === */
/* debug output folder */
#ifndef DEBUG_FOLDER
#define DEBUG_FOLDER "debugout-mmore"
#endif
/* root of project directory */
#ifndef PROJECT_LOC
#define PROJECT_LOC ./ project_loc /
#endif
/* location of binary */
#ifndef INSTALL_LOC
#define INSTALL_LOC ./ install_loc /
#endif
/* location of HMMER binary folder */
#ifndef HMMER_BIN_LOC
#define HMMER_BIN_LOC ./ hmmer_bin /
#endif
/* location of MMSEQS binary folder */
#ifndef MMSEQS_BIN_LOC
#define MMSEQS_BIN_LOC / mmseqs_bin /
#endif
/* location of MMSEQS binary folder */
#ifndef MMORE_BIN_LOC
#define MMORE_BIN_LOC / mmore_bin /
#endif

/* === SPECIAL BUILD OPTIONS === */
/* Debug - Extra debugging output */
#ifndef DEBUG
#define DEBUG FALSE
#endif
/* Memory Checks - Verifies all memory is clean before/after algorithms */
#ifndef MEMCHECK
#define MEMCHECK FALSE
#endif
/* Visualizer - Adds output for visualizing data */
#ifndef VIZ
#define VIZ DEBUG
#endif
/* Safe - Edgebound checking */
#ifndef SAFE
#define SAFE FALSE
#endif
