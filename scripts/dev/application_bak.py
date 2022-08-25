#!/usr/bin/env python
###############################################################################
#  - FILE: 	application.py
#  - DESC:  Interface for using MMORESEQS search.
###############################################################################

from email.policy import default
import os
import sys
from tokenize import Double, Double3
import click
import json
import typing
import functools

# Include mmoreseqs module path
py_mmoreseqs_module_path = os.getcwd() + "/../build"
sys.path.append(py_mmoreseqs_module_path)

# Try to import mmoreseqs
try:
    import py_mmoreseqs
    try:
        ans = py_mmoreseqs.add(4, 3)
        print("ans =", ans)
    except:
        print("ERROR: py_mmoreseqs package is not compiled properly.")
except:
    print("ERROR: py_mmoreseqs package not found")

opt_types = {
    # Debug options
    "--dbg": bool,
    "--dbg-viz": bool,
    # Pipeline options
    "--verbose": int,
    "--num-threads": int,
    "--search-type": str,
    # Task options
    "--run-mmseqs": bool,
    "--run-mmseqs-pref": bool,
    "--run-mmseqs-align": bool,
    "--run-convert": bool,
    "--run-mmore": bool,
    "--use-pvals": bool,
    # Input programs
    "--program-mmseqs": str,
    "--program-hmmer": str,
    "--program-mmoreseqs": str,
    "--script-dir": str,
    # Input files
    "--index": str,
    "--local-tools": str,
    "--guess-ftype": bool,
    "--mmoreseqs-ftype": str,
    "--mmoreseqs-main-ftype": str,
    "--tmp": str,
    "--prep": str,
    "--mmseqs-m8": str,
    # Input data
    "--dbsizes": (Double, Double),
    # MMORE parameters
    "--alpha": Double,
    "--beta": Double,
    "--gamma": Double,
    "--hard-limit": Double,
    # MMORE options
    "--run-prep": bool,
    "--run-bias": bool,
    "--run-full": bool,
    "--run-domains": bool,
    "--run-vit-mmore": bool,
    "--run-mmseqsaln": bool,
    "--run-vitaln": bool,
    "--run-vit": bool,
    "--run-postaln": bool,
    "--run-post": bool,
    # MMORE filters
    "--run-filter": bool,
    "--run-vit-filter": bool,
    "--run-cld-filter": bool,
    "--run-fwd-filter": bool,
    "--vit-filter": Double,
    "--cld-filter": Double,
    "--fwd-filter": Double,
    "--eval": Double,
    # MMSEQS parameters
    "--mmseqs-split": int,
    "--mmseqs-kmer": int,
    "--mmseqs-kscore": int,
    "--mmseqs-sens": Double,
    "--mmseqs-ungapped-vit": Double,
    "--mmseqs-eval": Double,
    "--mmseqs-pval": Double,
    "--mmseqs-hits-per-search": int,
    "--mmseqs-altalis": int,
    # MMSEQs data
    "--mmseqs-times": Double,
    "--mmseqs-dbsizes": int,
    # Search range output
    "--range": (int, int),
    "--search-mode": int,
    # Output
    "--mmseqs-m8out": str,
    "--stderr": str,
    "--stdout": str,
    "--allout": str,
    "--domtblout": str,
    "--m8out": str,
    "--myout": str,
    "--mydomtblout": str,
    "--mytimeout": str,
    "--mythreshout": str,
    "--output": str,
    "--customout": str,
    "--debugout": str
}

opt_defaults = {
    # Debug options
    "--dbg": False,
    "--dbg-viz": False,
    # Pipeline options
    "--verbose": 3,
    "--num-threads": 1,
    "--search-type": "P2S",
    # Task options
    "--run-mmseqs": True,
    "--run-mmseqs-pref": True,
    "--run-mmseqs-align": True,
    "--run-convert": True,
    "--run-mmore": True,
    "--use-pvals": bool,
    # Input programs
    "--program-mmseqs": None,
    "--program-hmmer": None,
    "--program-mmoreseqs": None,
    "--script-dir": None,
    # Input files
    "--index": None,
    "--local-tools": None,
    "--guess-ftype": None,
    "--mmoreseqs-ftype": str,
    "--mmoreseqs-main-ftype": str,
    "--tmp": str,
    "--prep": True,
    "--mmseqs-m8": True,
    # Input data
    "--dbsizes": True,
    # MMORE parameters
    "--alpha": Double,
    "--beta": Double,
    "--gamma": Double,
    "--hard-limit": Double,
    # MMORE options
    "--run-prep": bool,
    "--run-bias": bool,
    "--run-full": bool,
    "--run-domains": bool,
    "--run-vit-mmore": bool,
    "--run-mmseqsaln": bool,
    "--run-vitaln": bool,
    "--run-vit": bool,
    "--run-postaln": bool,
    "--run-post": bool,
    # MMORE filters
    "--run-filter": bool,
    "--run-vit-filter": bool,
    "--run-cld-filter": bool,
    "--run-fwd-filter": bool,
    "--vit-filter": Double,
    "--cld-filter": Double,
    "--fwd-filter": Double,
    "--eval": Double,
    # MMSEQS parameters
    "--mmseqs-split": int,
    "--mmseqs-kmer": int,
    "--mmseqs-kscore": int,
    "--mmseqs-sens": Double,
    "--mmseqs-ungapped-vit": Double,
    "--mmseqs-eval": Double,
    "--mmseqs-pval": Double,
    "--mmseqs-hits-per-search": int,
    "--mmseqs-altalis": int,
    # MMSEQs data
    "--mmseqs-times": Double,
    "--mmseqs-dbsizes": int,
    # Search range output
    "--range": (int, int),
    "--search-mode": int,
    # Interrim output
    "--mmseqs-m8out": None,
    # Output
    "--stderr": None,
    "--stdout": None,
    "--allout": None,
    "--domtblout": None,
    "--m8out": None,
    "--myout": None,
    "--mydomtblout": None,
    "--mytimeout": None,
    "--mythreshout": None,
    "--output": None,
    "--customout": None,
    "--debugout": None
}

opt_help = {
    # Debug options
    "--dbg": "Outputs debugger output.",
    "--dbg-viz": "Outputs visualization debugger output.",
    # Pipeline options
    "--verbose": "Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info",
    "--num-threads": "Number of parallel threads.",
    "--search-type": "Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence",
    # Task options
    "--run-mmseqs": "Run MMSEQS stage of MMORESEQS pipeline.",
    "--run-mmseqs-pref": "Run MMSEQS prefilter for MMSEQS stage of pipeline.",
    "--run-mmseqs-align": "Run Viterbi alignment for MMSEQS stage of pipeline.",
    "--run-convert": "Convert results from MMSEQS to a format that can be interpreted by MMORE.",
    "--run-mmore": "Run MMORE stage of the MMORESEQS pipeline.",
    "--use-pvals": "Use P-values (as opposed to E-values) for reporting and scoring thresholds.",
    # Input programs
    "--program-mmseqs": "Specify location of MMSEQS program.",
    "--program-hmmer": "Specify location of HMMER program.",
    "--program-mmoreseqs": "Specify location of MMORESEQS program.",
    "--script-dir": "Specify location of MMORESEQS scripts directory.",
    # Input files
    "--index": "Index files for target and query files.",
    "--local-tools": "Specify whether to use locally tools or globally installed tools.",
    "--guess-ftype": "Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.",
    "--mmoreseqs-ftype": "Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.",
    "--mmoreseqs-main-ftype": "Specify filetype of [0] target and [1] query input.",
    "--temp": "Specify temporary working folder location.",
    "--mmseqs-m8": "Specify location of input MMSEQS .m8 results file.",
    # Input data
    "--dbsizes": "Specify size of [0] target and [1] query databases.",
    # MMORE parameters
    "--alpha": "Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.",
    "--beta": "Set MMORE beta parameter. Determines X-drop below global max score accepted.",
    "--gamma": "Set MMORE gamma parameter. Determines X-drop below global max score before ending search.",
    "--hard-limit": "Set MMORE hard limit parameter. Determines lowest permissable score before ending search.",
    # MMORE options
    "--run-prep": "Run file preparation for the MMORE stage of pipeline",
    "--run-bias": "Use null bias during the MMORE stage of pipeline",
    "--run-full": "Run full quadratic search during the MMORE stage of pipeline. Overrides alpha/beta/gamma parameters.",
    "--run-domains": "Run over all domains found by MMSEQS for all target-query pairs during MMORE stage of pipeline. Alternatively, only uses the searches the top domain.",
    "--run-vit-mmore": "Run full viterbi during the MMORE stage of pipeline.",
    "--run-vit": "Run viterbi during the MMORE stage of pipeline.",
    "--run-post": "Run posterior during the MMORE stage of pipeline.",
    "--range": "Specify [0] start and [1] stop range of targets to search against the query database.",
    "--search-mode": "Specify mode of search. By default, searches entire query db against entire target db.",
    # MMORE filters
    "--run-filter": "Run all filters during the MMORE stage of pipeline .",
    "--run-vit-filter": "Run viterbi filter during the MMORE stage of pipeline.",
    "--run-cld-filter": "Run cloud filter during the MMORE stage of pipeline.",
    "--run-fwd-filter": "Run forward filter during the MMORE stage of pipeline.",
    "--vit-filter": "Set viterbi filter threshold cutoff score.",
    "--cld-filter": "Set cloud filter threshold cutoff score.",
    "--fwd-filter": "Set forward filter threshold cutoff score.",
    "--eval": "Set E-value filter threshold cutoff score for reporting.",
    # MMSEQS parameters
    "--mmseqs-split": "Set MMSEQS split size to number of files to break database into.",
    "--mmseqs-kmer": "Set MMSEQS kmer search length.",
    "--mmseqs-kscore": "Set MMSEQS k-score filter threshold score.",
    "--mmseqs-sens": "Set MMSEQS sensitivity filter threshold scpre. Overrides k-score filter.",
    "--mmseqs-ungapped-vit": "Set MMSEQS ungapped viterbi filter threshold.",
    "--mmseqs-eval": "Set MMSEQS E-value filter threshold.",
    "--mmseqs-pval": "Set MMSEQS P-value filter threshold. Overrides E-value.",
    "--mmseqs-hits-per-search": "Sets max number of separate alignments to pass filter for each target/query pair.",
    "--mmseqs-altalis": "Report alternate alignments during the MMSEQS stage of pipeline.",
    # MMSEQs data
    "--mmseqs-times": "Specify MMSEQS runtime.",
    "--mmseqs-dbsizes": "Specify MMSEQS [0] target and [1] query database sizes.",
    # Alignment reporting.
    "--run-mmseqsaln": "Recover and report alignments during the MMSEQS stage of pipeline.",
    "--run-vitaln": "Recover the report the viterbi alignments during the MMORE stage of pipeline.",
    "--run-postaln": "Recover and report the posterior alignments during the MMORE stage of pipeline.",
    # Output
    "--mmseqs-m8out": "Specify path to output MMSEQS results to .m8 format tsv file.",
    "--stderr": "Specify path to redirect standard error output.",
    "--stdout": "Specify path to redirect standard output",
    "--allout": "Specify path to output MMORESEQS results to exhaustive tsv file.",
    "--domtblout": "Specify path to output MMORESEQS results to HMMER-style domain table tsv file.",
    "--m8out": "Specify path to output MMORESEQS results to .m8 format tsv file.",
    "--myout": "Specify path to output MMORESEQS results to custom tsv file.",
    "--mydomtblout": "Specify path to output MMORESEQS results to custom domain table tsv file.",
    "--mytimeout": "Specify path to output MMORESEQS runtimes to tsv file.",
    "--mythreshout": "Specify path to output MMORESEQS threshold scores to tsv file.",
    "--customout": "Specify path to output MMORESEQS custom output to file.",
    "--debugout": "Specify path to output debugging data."
}

opt_hidden = {
    # Debug options
    "--dbg": True,
    "--dbg-viz": True,
    # Pipeline options
    "--verbose": False,
    "--num-threads": False,
    "--search-type": False,
    # Task options
    "--run-mmseqs": True,
    "--run-mmseqs-pref": True,
    "--run-mmseqs-align": True,
    "--run-convert": True,
    "--run-mmore": True,
    "--use-pvals": True,
    # Input programs
    "--program-mmseqs": True,
    "--program-hmmer": True,
    "--program-mmoreseqs": True,
    "--script-dir": True,
    # Input files
    "--index": True,
    "--local-tools": False,
    "--guess-ftype": False,
    "--mmoreseqs-ftype": False,
    "--mmoreseqs-main-ftype": False,
    "--tmp": False,
    "--prep": True,
    "--mmseqs-m8": True,
    # Input data
    "--dbsizes": True,
    # MMORE parameters
    "--alpha": True,
    "--beta": True,
    "--gamma": True,
    "--hard-limit": True,
    # MMORE options
    "--run-prep": True,
    "--run-bias": True,
    "--run-full": True,
    "--run-domains": True,
    "--run-vit-mmore": True,
    "--run-mmseqsaln": True,
    "--run-vitaln": True,
    "--run-vit": True,
    "--run-postaln": True,
    "--run-post": True,
    # MMORE filters
    "--run-filter": True,
    "--run-vit-filter": True,
    "--run-cld-filter": True,
    "--run-fwd-filter": True,
    "--vit-filter": True,
    "--cld-filter": True,
    "--fwd-filter": True,
    "--eval": True,
    # MMSEQS parameters
    "--mmseqs-split": True,
    "--mmseqs-kmer": True,
    "--mmseqs-kscore": True,
    "--mmseqs-sens": True,
    "--mmseqs-ungapped-vit": True,
    "--mmseqs-eval": True,
    "--mmseqs-pval": True,
    "--mmseqs-hits-per-search": True,
    "--mmseqs-altalis": True,
    # MMSEQs data
    "--mmseqs-times": True,
    "--mmseqs-dbsizes": True,
    # Search range output
    "--range": True,
    "--search-mode": True,
    # Output
    "--mmseqs-m8out": True,
    "--stderr": True,
    "--stdout": True,
    "--allout": True,
    "--domtblout": True,
    "--m8out": True,
    "--myout": True,
    "--mydomtblout": True,
    "--mytimeout": True,
    "--mythreshout": True,
    "--output": True,
    "--customout": True,
    "--debugout": True
}

opts = opt_defaults


def add_options(opt_names, func):
    """
    Add all options in list as decorations to command function.
    """
    for i in range(len(opt_names)):
        func = add_option(opt_names[i], func)
    return func


def add_option(opt_name, func):
    """
    Adds option with given name as decorator to command function.
    """
    # @click.option(opt_name, is_flag=(opt_types[opt_name] == bool), default=opt_defaults[opt_name], help=opt_help[opt_name], hidden=opt_hidden[opt_name])
    @click.option(opt_name)
    def set_option(opt_name):
        print("set_option: ", opt_name)
        return

    def wrapper_common_options(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper_common_options


def common_options(func):
    """
    Common arguments for multiple commands.
    """
    @click.option("--dbg")
    # @click.option("--dbg-viz", is_flag=True, default=False)(func)
    # @click.option("--verbose", nargs=1, default=0)(func)
    # @click.option("--num-threads")
    @click.option("--search-type")
    @functools.wraps(func)
    def wrapper_common_options(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper_common_options


@click.group()
def cli():
    """
    MMOREseqs command suite.
    """


@cli.command()
def version():
    """
    Display the current software version.
    """
    version = "1.0."
    print("VERSION: ", version)


@cli.command()
@click.argument("fasta_file", type=str)
@click.argument("hmm_file", type=str)
@common_options
def prep(fasta_file: str, hmm_file: str):
    """
    Preps input files for MMOREseqs search.
    """
    print("MMORESEQS PREP")
    print("FASTA_FILE: ", fasta_file)
    print("HMM_FILE: ", hmm_file)


@cli.command()
@click.argument("fasta_file", type=str)
@click.argument("hmm_file", type=str)
@click.argument("temp_file", type=str)
@common_options
def mmseqs_search(fasta_file: str, hmm_file: str, temp_file: str):
    """
    Perform MMSEQS stage of MMOREseqs search.
    """
    print("MMORESEQS MMSEQS_SEARCH")
    print("FASTA_FILE: ", fasta_file)
    print("HMM_FILE: ", hmm_file)
    print("TEMP_FILE: ", temp_file)


@cli.command()
@click.argument("fasta_file")
@click.argument("hmm_file")
@click.argument("temp_file")
@common_options
def mmore_search(fasta_file: str, hmm_file: str, temp_file: str):
    """
    Perform MMORE stage of MMOREseqs search.
    """
    print("MMORESEQS MMORE_SEARCH")
    print("FASTA_FILE: ", fasta_file)
    print("HMM_FILE: ", hmm_file)
    print("TEMP_FILE: ", temp_file)


@cli.command()
@click.argument("fasta_file")
@click.argument("hmm_file")
@click.argument("temp_file")
@common_options
def easy_search(fasta_file: str, hmm_file: str, temp_file: str):
    """
    Perform complete MMOREseqs search.
    """
    print("MMORESEQS MMORE_SEARCH")
    print("FASTA_FILE: ", fasta_file)
    print("HMM_FILE: ", hmm_file)
    print("TEMP_FILE: ", temp_file)


if __name__ == '__main__':
    cli()
