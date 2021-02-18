AUTOMAKE_OPTIONS = foreign subdir-objects

bin_PROGRAMS = mmore

mmore_CFLAGS = \
    -I./lib/easel/ \
    -L./lib/easel/ 

mmore_LDADD = \
    -leasel \
	 -lm 

mmore_SOURCES = \
    src/application.c \
    src/algs_linear/bound_fwdbck_linear.c \
    src/algs_linear/cloud_search_linear.c \
    src/algs_linear/fwdback_linear.c \
    src/algs_linear/merge_reorient_linear.c \
    src/algs_linear/pruning_linear.c \
    src/algs_linear/traceback_linear.c \
    src/algs_linear/viterbi_linear.c \
    src/algs_quad/bound_fwdbck_quad.c \
    src/algs_quad/bound_viterbi_quad.c \
    src/algs_quad/cloud_search_quad.c \
    src/algs_quad/fwdback_quad.c \
    src/algs_quad/merge_reorient_quad.c \
    src/algs_quad/posterior_quad.c \
    src/algs_quad/pruning_quad.c \
    src/algs_quad/ranged_fwdback_quad.c \
    src/algs_quad/traceback_quad.c \
    src/algs_quad/viterbi_quad.c \
    src/algs_sparse/bound_fwdbck_sparse.c \
    src/algs_sparse/bound_viterbi_sparse.c \
    src/algs_sparse/viterbi_traceback_sparse.c \
    src/algs_sparse/posterior_sparse.c \
    src/algs_sparse/posterior_null2_sparse.c \
    src/algs_sparse/posterior_optacc_sparse.c \
    src/algs_sparse/posterior_traceback_sparse.c \
    src/objects/alignment.c \
    src/objects/args.c \
    src/objects/clock.c \
    src/objects/debugger.c \
    src/objects/domain_def.c \
    src/objects/domain_x.c \
    src/objects/dp_matrix.c \
    src/objects/edgebound_rows.c \
    src/objects/edgebound.c \
    src/objects/f_index.c \
    src/objects/hmm_bg.c \
    src/objects/hmm_profile.c \
    src/objects/mystring.c \
    src/objects/results.c \
    src/objects/score_matrix.c \
    src/objects/sequence.c \
    src/objects/structs.c \
    src/objects/worker.c \
    src/objects/worker_thread.c \
    src/objects/io/reader.c \
    src/objects/basic/bound.c \
    src/objects/basic/char.c \
    src/objects/basic/double.c \
    src/objects/basic/float.c \
    src/objects/basic/int.c \
    src/objects/basic/range.c \
    src/objects/basic/str.c \
    src/objects/basic/template.c \
    src/objects/basic/trace.c \
    src/objects/basic/x_string.c \
    src/objects/matrix/matrix_2d.c \
    src/objects/matrix/matrix_3d.c \
    src/objects/matrix/matrix_3d_sparse.c \
    src/objects/vectors/vector_bound.c \
    src/objects/vectors/vector_char.c \
    src/objects/vectors/vector_double.c \
    src/objects/vectors/vector_float.c \
    src/objects/vectors/vector_int.c \
    src/objects/vectors/vector_range.c \
    src/objects/vectors/vector_template.c \
    src/objects/vectors/vector_trace.c \
    src/parsers/arg_parser.c \
    src/parsers/hmm_parser.c \
    src/parsers/index_parser.c \
    src/parsers/m8_parser.c \
    src/parsers/seq_parser.c \
    src/parsers/seq_to_profile.c \
    src/pipelines/pipeline_hmmbuild.c \
    src/pipelines/pipeline_index.c \
    src/pipelines/pipeline_itest.c \
    src/pipelines/pipeline_main.c \
    src/pipelines/pipeline_mmseqs.c \
    src/pipelines/pipeline_mmseqs_plus.c \
    src/pipelines/pipeline_null.c \
    src/pipelines/pipeline_time.c \
    src/pipelines/pipeline_utest.c \
    src/pipelines/pipeline_visualize.c \
    src/pipelines/report.c \
    src/pipelines/work.c \
    src/utilities/error_handler.c \
    src/utilities/dp_matrix_util.c \
    src/utilities/misc.c \
    src/utilities/rng.c \
    src/utilities/testing.c

mmore_HEADERS = \
\