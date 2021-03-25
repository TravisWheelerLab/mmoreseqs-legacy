cmake_minimum_required(VERSION 3.10)
include(ExternalProject)

project(mmoreseqs VERSION 1.0.0)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)

# -----------
# Executables
# -----------

add_executable(
        mmoreseqs
        src/application.c

        src/algs_linear/bound_fwdbck_linear.c
        src/algs_linear/cloud_search_linear.c
        src/algs_linear/pruning_linear.c
        src/algs_linear/viterbi_linear.c
        src/algs_linear/viterbi_traceback_linear.c
        src/algs_quad/bound_fwdbck_quad.c
        src/algs_quad/cloud_search_quad.c
        src/algs_quad/fwdback_quad.c
        src/algs_quad/posterior_quad.c
        src/algs_quad/pruning_quad.c
        src/algs_quad/viterbi_quad.c
        src/algs_quad/viterbi_traceback_quad.c
        src/algs_sparse/bound_fwdbck_sparse.c
        src/algs_sparse/bound_viterbi_sparse.c
        src/algs_sparse/posterior_null2_sparse.c
        src/algs_sparse/posterior_sparse.c
        src/algs_sparse/posterior_traceback_sparse.c
        src/algs_sparse/viterbi_traceback_sparse.c
        src/objects/alignment.c
        src/objects/arg_opts.c
        src/objects/args.c
        src/objects/basic/bound.c
        src/objects/basic/char.c
        src/objects/basic/double.c
        src/objects/basic/float.c
        src/objects/basic/gen.c
        src/objects/basic/gen_ext.c
        src/objects/basic/int.c
        src/objects/basic/ptr.c
        src/objects/basic/range.c
        src/objects/basic/str.c
        src/objects/basic/str_ext.c
        src/objects/basic/trace.c
        src/objects/clock.c
        src/objects/commandline.c
        src/objects/domain_def.c
        src/objects/f_index.c
        src/objects/hmm_bg.c
        src/objects/hmm_profile.c
        src/objects/io/buffer.c
        src/objects/io/buffer_read.c
        src/objects/io/filer.c
        src/objects/io/reader.c
        src/objects/m8_results.c
        src/objects/matrix/matrix_2d.c
        src/objects/matrix/matrix_3d.c
        src/objects/matrix_sparse/edgebound.c
        src/objects/matrix_sparse/edgebound_merge_reorient.c
        src/objects/matrix_sparse/edgebound_rows.c
        src/objects/matrix_sparse/matrix_3d_sparse.c
        src/objects/matrix_sparse/matrix_3d_sparse_build.c
        src/objects/mystring.c
        src/objects/results.c
        src/objects/score_matrix.c
        src/objects/scriptrunner.c
        src/objects/sequence.c
        src/objects/structs.c
        src/objects/vectors/vector_bound.c
        src/objects/vectors/vector_char.c
        src/objects/vectors/vector_double.c
        src/objects/vectors/vector_float.c
        src/objects/vectors/vector_int.c
        src/objects/vectors/vector_ptr.c
        src/objects/vectors/vector_range.c
        src/objects/vectors/vector_str.c
        src/objects/vectors/vector_trace.c
        src/objects/worker.c
        src/parsers/arg_parser.c
        src/parsers/hmm_parser.c
        src/parsers/index_parser.c
        src/parsers/m8_parser.c
        src/parsers/seq_parser.c
        src/parsers/seq_to_profile.c
        src/parsers/seq_to_profile.c
        src/pipelines/pipeline_generic.c
        src/pipelines/pipeline_hmmbuild.c
        src/pipelines/pipeline_index.c
        src/pipelines/pipeline_interactive.c
        src/pipelines/pipeline_mmoreseqs_easysearch.c
        src/pipelines/pipeline_mmoreseqs_mmore.c
        src/pipelines/pipeline_mmoreseqs_mmseqs.c
        src/pipelines/pipeline_mmoreseqs_prep.c
        src/pipelines/pipeline_mmoreseqs_prepsearch.c
        src/pipelines/pipeline_mmoreseqs_search.c
        src/pipelines/pipeline_null.c
        src/reporting/domtblout.c
        src/reporting/m8out.c
        src/reporting/mainout.c
        src/reporting/mydomout.c
        src/reporting/mydomout.c
        src/reporting/myout.c
        src/reporting/mythreshout.c
        src/reporting/mytimeout.c
        src/reporting/report_util.c
        src/utilities/error_handler.c
        src/utilities/mymath.c
        src/utilities/rng.c
        src/utilities/statistics.c
        src/utilities/system_io.c
        src/work/work_cloud_fwdback.c
        src/work/work_cloud_merge.c
        src/work/work_cloud_search.c
        src/work/work_etc.c
        src/work/work_index.c
        src/work/work_loader.c
        src/work/work_loop.c
        src/work/work_maintenance.c
        src/work/work_optacc.c
        src/work/work_posterior.c
        src/work/work_posterior_bydom.c
        src/work/work_report.c
        src/work/work_scoring.c
        src/work/work_scripting.c
        src/work/work_sparse_mx.c
        src/work/work_threshold.c
        src/work/work_viterbi.c
)

# ------------
# Dependencies
# ------------

ExternalProject_Add(
        easel
        URL https://github.com/EddyRivasLab/easel/archive/Bio-Easel-0.14.tar.gz
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/easel
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoconf
        COMMAND ./configure --prefix=${SOURCE_DIR}
        BUILD_COMMAND "make"
        INSTALL_COMMAND ""
        TEST_COMMAND ""
)
ExternalProject_Get_property(easel SOURCE_DIR)
target_include_directories(mmoreseqs PRIVATE ${SOURCE_DIR})
target_link_libraries(mmoreseqs ${SOURCE_DIR}/libeasel.a)
add_dependencies(mmoreseqs easel)

# Add math dependency
target_link_libraries(mmoreseqs m)

