cmake_minimum_required(VERSION 3.10)

# disallow building in same directory as source.
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(FATAL_ERROR "ERROR: Cmake build directory cannot be the same as root source directory. Hint: Can build with `make`.")
endif()

# set(CMAKE_DISABLE_SOURCE_CHANGES ON)
# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
project(mmoreseqs VERSION "1.0.0")

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

include(ExternalProject)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -----------
# Input Files
# -----------
set(SOURCE_FOLDER src)
set(SOURCES
  src/objects/sequence.c
  src/objects/mystring.c
  src/objects/domain_x.c
  src/objects/arg_opts.c
  src/objects/f_index.c
  src/objects/m8_results.c
  src/objects/worker_thread.c
  src/objects/score_matrix.c
  src/objects/alignment.c
  src/objects/commandline.c
  src/objects/matrix_sparse/edgebound.c
  src/objects/matrix_sparse/matrix_3d_sparse_build.c
  src/objects/matrix_sparse/edgebound_merge_reorient.c
  src/objects/matrix_sparse/matrix_3d_sparse.c
  src/objects/matrix_sparse/edgebound_rows.c
  src/objects/hmm_profile.c
  src/objects/scriptrunner.c
  src/objects/worker.c
  src/objects/io/writer.c
  src/objects/io/buffer_write.c
  src/objects/io/reader.c
  src/objects/io/buffer_read.c
  src/objects/io/buffer.c
  src/objects/io/filer.c
  src/objects/basic/bound.c
  src/objects/basic/ptr.c
  src/objects/basic/xstr.c
  src/objects/basic/gen_ext.c
  src/objects/basic/char.c
  src/objects/basic/bool.c
  src/objects/basic/double.c
  src/objects/basic/str.c
  src/objects/basic/gen.c
  src/objects/basic/template.c
  src/objects/basic/float.c
  src/objects/basic/range.c
  src/objects/basic/str_ext.c
  src/objects/basic/int.c
  src/objects/basic/trace.c
  src/objects/debugger.c
  src/objects/matrix/matrix_2d.c
  src/objects/matrix/matrix_3d.c
  src/objects/domain_def.c
  src/objects/results.c
  src/objects/map/xxx_umap.c
  src/objects/args.c
  src/objects/dp_matrix.c
  src/objects/vectors/vector_float.c
  src/objects/vectors/vector_bound.c
  src/objects/vectors/vector_double.c
  src/objects/vectors/vector_trace.c
  src/objects/vectors/vector_ptr.c
  src/objects/vectors/vector_int.c
  src/objects/vectors/vector_char.c
  src/objects/vectors/vector_template.c
  src/objects/vectors/vector_str.c
  src/objects/vectors/vector_range.c
  src/objects/clock.c
  src/objects/structs.c
  src/objects/hmm_bg.c
  src/reporting/hmmerout.c
  src/reporting/domtblout.c
  src/reporting/mydomout.c
  src/reporting/m8out.c
  src/reporting/mythreshout.c
  src/reporting/mainout.c
  src/reporting/report_util.c
  src/reporting/mytimeout.c
  src/reporting/myout.c
  src/work/work_optacc.c
  src/work/work_viterbi.c
  src/work/work_loader.c
  src/work/work_report.c
  src/work/work_posterior_bydom.c
  src/work/work_cloud_search.c
  src/work/work_posterior.c
  src/work/work_sparse_mx.c
  src/work/work_loop.c
  src/work/work_cloud_fwdback.c
  src/work/work_maintenance.c
  src/work/work_etc.c
  src/work/work_index.c
  src/work/work_scoring.c
  src/work/work_fwdback.c
  src/work/work_scripting.c
  src/work/work_threshold.c
  src/work/work_cloud_merge.c
  src/application.c
  src/parsers/seq_parser.c
  src/parsers/arg_parser.c
  src/parsers/m8_parser.c
  src/parsers/seq_to_profile.c
  src/parsers/hmm_parser.c
  src/parsers/index_parser.c
  src/algs_naive/bound_fwdbck_naive.c
  src/algs_vectorized/fwdback_vectorized.c
  src/algs_quad/pruning_quad.c
  src/algs_quad/bound_posterior_quad.c
  src/algs_quad/viterbi_traceback_quad.c
  src/algs_quad/bound_viterbi_quad.c
  src/algs_quad/bound_fwdbck_quad.c
  src/algs_quad/merge_reorient_quad.c
  src/algs_quad/viterbi_quad.c
  src/algs_quad/fwdback_quad.c
  src/algs_quad/posterior_quad.c
  src/algs_quad/cloud_search_quad.c
  src/pipelines/workflow.c
  src/pipelines/pipeline_mmoreseqs_mmore.c
  src/pipelines/pipeline_mmoreseqs_prepsearch.c
  src/pipelines/pipeline_mmoreseqs_easysearch.c
  src/pipelines/pipeline_null.c
  src/pipelines/pipeline_mmoreseqs_mmseqs.c
  src/pipelines/pipeline_mmoreseqs_prep.c
  src/pipelines/pipeline_index.c
  src/pipelines/pipeline_utest.c
  src/pipelines/pipeline_mmoreseqs_search.c
  src/algs_sparse/bound_fwdbck_sparse_test.c
  src/algs_sparse/bound_fwdbck_sparse.c
  src/algs_sparse/viterbi_traceback_sparse.c
  src/algs_sparse/bound_viterbi_sparse.c
  src/algs_sparse/posterior_traceback_sparse.c
  src/algs_sparse/posterior_null2_sparse.c
  src/algs_sparse/posterior_sparse.c
  src/algs_linear/bound_fwdbck_linear.c
  src/algs_linear/fwdback_linear.c
  src/algs_linear/viterbi_traceback_linear.c
  src/algs_linear/pruning_linear.c
  src/algs_linear/viterbi_linear.c
  src/algs_linear/cloud_search_linear.c
  src/utilities/testing.c
  src/utilities/visualizer.c
  src/utilities/error_handler.c
  src/utilities/system_io.c
  src/utilities/rng.c
  src/utilities/mymath.c
  src/utilities/dp_matrix_util.c
  src/utilities/statistics.c
)

set(SCRIPT_FOLDER scripts)
set(SCRIPTS
  scripts/workflows/mmoreseqs-prep.sh
  scripts/workflows/mmoreseqs-prepsearch.sh
  scripts/workflows/mmoreseqs-mmseqssearch.sh
  scripts/workflows/mmoreseqs-mmoresearch.sh
  scripts/workflows/mmoreseqs-easysearch.sh
  scripts/workflows/helpers/translate_ids.sh
  scripts/workflows/helpers/translate_ids/make_dummy_align_mmdb.py
  scripts/workflows/helpers/mmore_functions.sh
  scripts/workflows/helpers/convert_fasta-to-hmm.sh
  scripts/workflows/mmoreseqs-search.sh
)

# ----------
# Copy Files
# ----------
# Python and Bash workflow scripts
add_subdirectory(${SCRIPT_FOLDER})
file(COPY ${SCRIPT_FOLDER} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# ------------
# Dependencies
# ------------

# Easel library
ExternalProject_Add(easel
  URL https://github.com/EddyRivasLab/easel/archive/Bio-Easel-0.14.tar.gz
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/easel
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND autoconf
  COMMAND ./configure --prefix=${SOURCE_DIR}
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)
ExternalProject_Get_property(easel SOURCE_DIR)
include_directories(BEFORE ${SOURCE_DIR})
link_libraries(${SOURCE_DIR}/libeasel.a)

# Math library
link_libraries(m)

# ---------
# Libraries
# ---------
# mmoreseqs library
add_library(
  mmoreseqs_lib
  SHARED
  ${SOURCES}
)
add_dependencies(mmoreseqs_lib easel)

# -----------
# Executables
# -----------

# mmoreseqs
add_executable(mmoreseqs
  src/mmoreseqs.c
  ${SOURCES}
)
add_dependencies(mmoreseqs easel)

# --------------------
# Build Configurations
# --------------------

# Options for setable macros
option(SET_BUILD_HASH "Set git commit hash of this build." "NO COMMIT HASH")
option(SET_HMMER_BIN "Set location of HMMER binary to use." "hmmer")
option(SET_MMSEQS_BIN "Set location of MMseqs binary to use." "mmseqs")
option(SET_DEBUG_FOLDER "Folder to send debug output files." "debug_out")
option(SET_DEBUG, "Set whether to output debug info." "0")

# Set macros
add_compile_definitions("BUILD_VERSION=\"${CMAKE_PROJECT_VERSION}\"")
add_compile_definitions("BUILD_TYPE=\"${CMAKE_BUILD_TYPE}\"")
add_compile_definitions("BUILD_HASH=\"${SET_BUILD_HASH}\"")
add_compile_definitions("BUILD_VERSION=\"${SET_BUILD_VERSION}\"")
add_compile_definitions("PROJECT_LOC=${CMAKE_CURRENT_BINARY_DIR}")
add_compile_definitions("SCRIPT_LOC=${CMAKE_CURRENT_BINARY_DIR}/scripts/")
add_compile_definitions("MMORE_BIN_LOC=${CMAKE_CURRENT_BINARY_DIR}/bin/")

# C compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG=1")
set(CMAKE_C_FLAGS_ASAN "-g -pg -fsanitize=address")

# CPP compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
