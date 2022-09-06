/*******************************************************************************
 *  - FILE:  cli_application.cpp
 *  - DESC:  Entry Point to Application. Parses commandline args and starts engine.
 *******************************************************************************/

#include "cli_application.hpp"

/* initialize worker and args object */
WORKER* worker;
ARGS* args;

int main(int argc, char* argv[]) {
  // worker = WORKER_Create();
  // worker = WORKER_Init(worker);
  // args = worker->args;

  CLI::App mmoreseqs_app{"MMOREseqs"};
  mmoreseqs_app.require_subcommand(1);
  Add_AllOptions(mmoreseqs_app);

  EasySearchCommand(mmoreseqs_app);
  PrepSearchCommand(mmoreseqs_app);
  SearchCommand(mmoreseqs_app);

  CLI11_PARSE(mmoreseqs_app, argc, argv);
  /* output arguments */
  // if (args->verbose_level >= VERBOSE_LOW) {
  //   ARGS_Dump(args, stdout);
  // }
  // if (args->pipeline_mode != PIPELINE_TEST) {
  //   PIPELINES[args->pipeline_mode].pipeline_main(worker);
  // } else {
  //   std::cout << "No valid subcommand selected." << std::endl;
  // }

  return EXIT_SUCCESS;
}

void EasySearchCommand(CLI::App& mmoreseqs_app) {
  CLI::App* easy_search_cmd = mmoreseqs_app.add_subcommand(
      "easy-search",
      "Runs complete MMOREseqs search pipeline on input files.");
  easy_search_cmd->add_option("query_msa")->required();
  easy_search_cmd->add_option("target_fasta")->required();
  easy_search_cmd->add_option("temp_dir")->required();
}

void PrepSearchCommand(CLI::App& mmoreseqs_app) {
  CLI::App* prep_cmd = mmoreseqs_app.add_subcommand(
      "prep",
      "Prepares input files for running MMOREseqs search and stores output in temporary directory.");
  prep_cmd->add_option("query_msa")->required();
  prep_cmd->add_option("target_fasta")->required();
  prep_cmd->add_option("prep_dir")->required();

  CLI::App* prep_search_cmd = mmoreseqs_app.add_subcommand(
      "prep-search",
      "Runs MMOREseqs search pipeline on prepared temporary directory.");
  prep_search_cmd->add_option("prep_dir")->required();
}

void SearchCommand(CLI::App& mmoreseqs_app) {
  CLI::App* search_cmd = mmoreseqs_app.add_subcommand(
      "search",
      "Runs MMOREseqs search pipeline on collection of query and target files. For internal and finetuned processes.");
  search_cmd->add_option("query_mmore_hmm")->required();
  search_cmd->add_option("target_mmore_fasta")->required();
  search_cmd->add_option("query_mmseqs_pmmdb")->required();
  search_cmd->add_option("query_mmseqs_smmdb")->required();
  search_cmd->add_option("target_mmseqs_smmdb")->required();

  CLI::App* mmore_search_cmd = mmoreseqs_app.add_subcommand(
      "mmore-search",
      "Runs MMORE phase of search pipeline on prepared temporary directory. For internal and finetuned processes.");
  mmore_search_cmd->add_option("query_mmore_hmm")->required();
  mmore_search_cmd->add_option("target_mmore_fasta")->required();

  CLI::App* mmseqs_search_cmd = mmoreseqs_app.add_subcommand(
      "mmseqs-search",
      "Runs MMSeqs search pipeline on prepared temporary directory. For internal and finetuned processes.");
  mmseqs_search_cmd->add_option("query_mmseqs_pmmdb")->required();
  mmseqs_search_cmd->add_option("query_mmseqs_smmdb")->required();
  mmseqs_search_cmd->add_option("target_mmseqs_smmdb")->required();
}

void Add_GeneralOptions(CLI::App& cmd) {
}

void Add_InternalOptions(CLI::App& cmd) {
}

void Add_DebugOptions(CLI::App& cmd) {
}

void Add_AllOptions(CLI::App& cmd) {
  CLI::Option_group* general_grp = cmd.add_option_group("General Options");
  general_grp->add_flag("--verbose/-v", "Amount of output. 0=quiet, 1=errors, 2=warning, 3=info");
  general_grp->add_flag("--num-threads/", "Number of parallel threads.");
  general_grp->add_flag("--eval/-E", "Set E-value filter threshold cutoff score for reporting.");
  general_grp->add_flag("--use-pvals/-P", "Use P-values (as opposed to default E-values) for reporting and scoring thresholds.");

  CLI::Option_group* pipeline_grp = cmd.add_option_group("Pipeline Options");
  pipeline_grp->add_flag("--run-prep", "Run file preparation for the MMORE stage of pipeline");
  pipeline_grp->add_flag("--run-mmseqs", "Run MMSEQS stage of MMORESEQS pipeline.");
  pipeline_grp->add_flag("--run-mmseqs-pref", "Run MMSEQS prefilter for MMSEQS stage of pipeline.");
  pipeline_grp->add_flag("--run-mmseqs-align", "Run Viterbi alignment for MMSEQS stage of pipeline.");
  pipeline_grp->add_flag("--run-convert", "Convert results from MMSEQS to a format that can be interpreted by MMORE.");
  pipeline_grp->add_flag("--run-mmore", "Run MMORE stage of the MMORESEQS pipeline.");

  CLI::Option_group* input_grp = cmd.add_option_group("Input Options");
  input_grp->add_flag("--temp/--prep", "Specify temporary working folder location.");
  input_grp->add_flag("--mmseqs-m8", "Specify location of input MMSEQS .m8 results file.");
  input_grp->add_flag("--index", "Specify location of [0] target and [1] query index files.");

  CLI::Option_group* reporting_grp = cmd.add_option_group("Reporting Options");
  reporting_grp->add_flag("--run-mmseqsaln", "Recover and report the viterbi alignments during the MMSEQS stage of pipeline.");
  reporting_grp->add_flag("--run-vitaln", "Recover the report the viterbi alignments during the MMORE stage of pipeline.");
  reporting_grp->add_flag("--run-postaln", "Recover and report the posterior alignments during the MMORE stage of pipeline.");

  CLI::Option_group* mmseqs_search_grp = cmd.add_option_group("MMSeqs Search Options");
  mmseqs_search_grp->add_flag("--mmseqs-kmer", "Sets the kmer length during the prefilter step of MMseqs stage.");
  mmseqs_search_grp->add_flag("--mmseqs-kscore", "Sets filter score for kmers during the prefilter step of MMseqs stage.");
  mmseqs_search_grp->add_flag("--mmseqs-hits-per-search", "Sets total number of hits per target allowed per search to be reported by prefilter step of MMseqs stage.");
  mmseqs_search_grp->add_flag("--mmseqs-ungapped-vit", "Sets filter score for ungapped viterbi step of MMseqs stage.");
  mmseqs_search_grp->add_flag("--mmseqs-eval", "Sets filter E-value of MMseqs for (gapped) viterbi step of MMseqs stage.");
  mmseqs_search_grp->add_flag("--mmseqs-pval", "Sets filter P-value of MMseqs for (gapped) viterbi step of MMseqs stage. Overrides E-value score.");
  mmseqs_search_grp->add_flag("--mmseqs-altalis", "Sets number of additional alternate alignments that can be found per target-query pair. Zero results in a single alignment.");
  mmseqs_search_grp->add_flag("--mmseqs-split", "When building MMseqs databases determines the number of independent files to store databases across. Larger number keeps file size smaller.");

  CLI::Option_group* mmore_search_grp = cmd.add_option_group("MMORE Search Options");
  mmore_search_grp->add_flag("--alpha", "Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.");
  mmore_search_grp->add_flag("--beta", "Set MMORE beta parameter. Determines X-drop below global max score accepted.");
  mmore_search_grp->add_flag("--gamma", "Set MMORE gamma parameter. Determines X-drop below global max score before ending search.");
  mmore_search_grp->add_flag("--hard-limit", "Set MMORE hard limit parameter. Determines lowest permissable score before ending search.");
  mmore_search_grp->add_flag("--vit-filter", "Set viterbi filter threshold cutoff score.");
  mmore_search_grp->add_flag("--cld-filter", "Set cloud filter threshold cutoff score.");
  mmore_search_grp->add_flag("--fwd-filter", "Set forward filter threshold cutoff score.");
  mmore_search_grp->add_flag("--range", "Specify [0] start and [1] stop range of targets to search against the query database.");
  mmore_search_grp->add_flag("--run-bias", "Use null bias during the MMORE stage of pipeline");
  mmore_search_grp->add_flag("--run-full", "Run full quadratic search during the MMORE stage of pipeline. Overrides alpha/beta/gamma parameters.");
  mmore_search_grp->add_flag("--run-domains", "Run search over all domains found by MMSEQS for all target-query pairs during MMORE stage of pipeline. Alternatively, only uses the searches the top domain.");
  mmore_search_grp->add_flag("--run-vit-mmore", "Run full viterbi during the MMORE stage of pipeline.");
  mmore_search_grp->add_flag("--run-vit", "Run viterbi during the MMORE stage of pipeline.");
  mmore_search_grp->add_flag("--run-post", "Run posterior during the MMORE stage of pipeline.");
  mmore_search_grp->add_flag("--run-filter", "Run all filters during the MMORE stage of pipeline.");
  mmore_search_grp->add_flag("--run-vit-filter", "Run viterbi filter during the MMORE stage of pipeline. Otherwise bypasses filter.");
  mmore_search_grp->add_flag("--run-cld-filter", "Run cloud filter during the MMORE stage of pipeline. Otherwise bypasses filter.");
  mmore_search_grp->add_flag("--run-fwd-filter", "Run forward filter during the MMORE stage of pipeline. Otherwise bypasses filter.");

  CLI::Option_group* forward_grp = cmd.add_option_group("Forwarding Options", "General users should never need to invoke these options. They store data which are passed to scripts for intermediary workflow operations.");
  forward_grp->add_flag("--search-type", "Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence");
  forward_grp->add_flag("--program-mmseqs", "Specify location of MMSEQS program.");
  forward_grp->add_flag("--program-hmmer", "Specify location of HMMER program.");
  forward_grp->add_flag("--program-mmoreseqs", "Specify location of MMORESEQS program.");
  forward_grp->add_flag("--script-dir", "Specify location of MMORESEQS scripts directory.");
  forward_grp->add_flag("--local-tools", "Specify whether to use locally tools or globally installed tools.");
  forward_grp->add_flag("--guess-ftype", "Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.");
  forward_grp->add_flag("--mmoreseqs-ftype", "Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.");
  forward_grp->add_flag("--mmoreseqs-main-ftype", "Specify filetype of [0] target and [1] query input.");
  forward_grp->add_flag("--dbsizes", "Specify size of [0] target and [1] query databases.");
  forward_grp->add_flag("--mmseqs-times", "Passthrough MMSEQS runtimes.");
  forward_grp->add_flag("--mmseqs-dbsizes", "Passthrough MMSEQS database sizes.");

  CLI::Option_group* dev_grp = cmd.add_option_group("Developer Options", "General users whould never need to invoke these options.  Used by developers for testing and debugging. Many require that project is built with `$ make build-debug`.");
  dev_grp->add_flag("--dbg", "Whether to output debugger information.");
  dev_grp->add_flag("--dbg-viz", "Whether to output debugger visualization information.");
  dev_grp->add_flag("--allout", "Outputs a exhaustive tsv file.");
  dev_grp->add_flag("--myout", "Outputs a custom formatted tsv file.");
  dev_grp->add_flag("--mytimeout", "Outputs a runtime tsv file.");
  dev_grp->add_flag("--mythreshout", "Outputs a threshold score tsv file.");
  dev_grp->add_flag("--myhmmerout", "Outputs a HMMER-style file.");
  dev_grp->add_flag("--debugout", "Output debugger information to file.");
}

void StoreFlagArgAsString(STR& flag_dest, std::string& flag_src) {
  STR_Set(flag_dest, flag_src.c_str());
}

void StoreFlagArgAsBool(bool& flag_dest, std::string& flag_src) {
  flag_dest = atoi(flag_src.c_str());
}

void StoreFlagArgAsInt(int& flag_dest, std::string& flag_src) {
  flag_dest = atoi(flag_src.c_str());
}

void StoreFlagArgAsDouble(int& flag_dest, std::string& flag_src) {
  flag_dest = atof(flag_src.c_str());
}
