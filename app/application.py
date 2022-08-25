#!/usr/bin/env python
###############################################################################
#  - FILE: 	application.py
#  - DESC:  Interface for using MMORESEQS search.
###############################################################################

# from email.policy import default
import os
import sys
from sre_constants import FAILURE
from tokenize import Double
import click
import json
import typing
import functools
import subprocess

# Try to import mmoreseqs
try:
    import mmoreseqs_pylib
    try:
        mmoreseqs_pylib.get_version()
    except Exception as err:
        print(err)
except Exception as err:
    print(err)

debug = False
primary_args = {
    "prep": ["target_filepath", "query_filepath", "temp_filepath"],
    "prep-search": ["temp_filepath"],
    "mmseqs-search": ["target_filepath", "query_filepath"],
    "mmore-search": ["target_filepath", "query_filepath"],
    "easy-search": ["target_filepath", "query_filepath", "temp_filepath"],
    "version": []
}


def dbg_print(*args, **kwargs):
    """
    Print if debugging. 
    """
    if debug:
        print(*args, **kwargs)

### CLICK COMMAND LINE ###


@click.group()
def cli():
    """
    MMOREseqs command suite.
    """
    pass


@cli.command(name="prep", help="Prepare files for performing MMORESEQS search via prep-search. Takes as arguments: [0] search target filepath, [1] search query filepath, and [2] temporary working directory path. Target filepath types: HMM. Query filepath types: HMM, FASTA.")
@click.argument("target_filepath", type=str)
@click.argument("query_filepath", type=str)
@click.argument("temp_filepath", type=str)
@click.option("--dbg", type=(bool), default=(None), help="Outputs debugger info.", hidden=True, nargs=1)
@click.option("--dbg-viz", type=(bool), default=(None), help="Outputs visualization debugger info.", hidden=True, nargs=1)
@click.option("--debugout", type=(bool), default=(None), help="Specify path to output debugger info.", hidden=True, nargs=1)
@click.option("--verbose", type=(int), default=(None), help="Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info", hidden=False, nargs=1)
@click.option("--num-threads", type=(int), default=(None), help="Number of parallel threads.", hidden=False, nargs=1)
@click.option("--search-type", type=(str), default=(None), help="Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence", hidden=False, nargs=1)
@click.option("--use-pvals", type=(bool), default=(None), help="Use P-values (as opposed to default E-values) for reporting and scoring thresholds.", hidden=False, nargs=1)
@click.option("--eval", type=(bool), default=(None), help="Set E-value filter threshold cutoff score for reporting.", hidden=False, nargs=1)
@click.option("--program-mmseqs", type=(str), default=(None), help="Specify location of MMSEQS program.", hidden=False, nargs=1)
@click.option("--program-hmmer", type=(str), default=(None), help="Specify location of HMMER program.", hidden=False, nargs=1)
@click.option("--program-mmoreseqs", type=(str), default=(None), help="Specify location of MMORESEQS program.", hidden=False, nargs=1)
@click.option("--script-dir", type=(str), default=(None), help="Specify location of MMORESEQS scripts directory.", hidden=False, nargs=1)
@click.option("--local-tools", type=(bool), default=(None), help="Specify whether to use locally tools or globally installed tools.", hidden=False, nargs=1)
@click.option("--guess-ftype", type=(bool), default=(None), help="Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=1)
@click.option("--mmoreseqs-ftype", type=(str, str, str, str), default=(None, None, None, None), help="Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=4)
def my_func(*args, **kwargs):
    prep(args, kwargs)
    pass


@cli.command(name="prep-search", help="Runs full MMORESEQS search on prepared file directory.  Takes as arguments: [0] temporary working directory path.")
@click.argument("temp_filepath", type=str)
@click.option("--dbg", type=(bool), default=(None), help="Outputs debugger info.", hidden=True, nargs=1)
@click.option("--dbg-viz", type=(bool), default=(None), help="Outputs visualization debugger info.", hidden=True, nargs=1)
@click.option("--debugout", type=(bool), default=(None), help="Specify path to output debugger info.", hidden=True, nargs=1)
@click.option("--verbose", type=(int), default=(None), help="Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info", hidden=False, nargs=1)
@click.option("--num-threads", type=(int), default=(None), help="Number of parallel threads.", hidden=False, nargs=1)
@click.option("--search-type", type=(str), default=(None), help="Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence", hidden=False, nargs=1)
@click.option("--run-mmseqs", type=(bool), default=(None), help="Run MMSEQS stage of MMORESEQS pipeline.", hidden=False, nargs=1)
@click.option("--run-mmseqs-pref", type=(bool), default=(None), help="Run MMSEQS prefilter for MMSEQS stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-mmseqs-align", type=(bool), default=(None), help="Run Viterbi alignment for MMSEQS stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-convert", type=(bool), default=(None), help="Convert results from MMSEQS to a format that can be interpreted by MMORE.", hidden=False, nargs=1)
@click.option("--run-mmore", type=(bool), default=(None), help="Run MMORE stage of the MMORESEQS pipeline.", hidden=False, nargs=1)
@click.option("--use-pvals", type=(bool), default=(None), help="Use P-values (as opposed to default E-values) for reporting and scoring thresholds.", hidden=False, nargs=1)
@click.option("--eval", type=(bool), default=(None), help="Set E-value filter threshold cutoff score for reporting.", hidden=False, nargs=1)
@click.option("--program-mmseqs", type=(str), default=(None), help="Specify location of MMSEQS program.", hidden=False, nargs=1)
@click.option("--program-hmmer", type=(str), default=(None), help="Specify location of HMMER program.", hidden=False, nargs=1)
@click.option("--program-mmoreseqs", type=(str), default=(None), help="Specify location of MMORESEQS program.", hidden=False, nargs=1)
@click.option("--script-dir", type=(str), default=(None), help="Specify location of MMORESEQS scripts directory.", hidden=False, nargs=1)
@click.option("--index", type=(str, str), default=(None, None), help="Specify location of [0] target and [1] query index files.", hidden=False, nargs=2)
@click.option("--local-tools", type=(bool), default=(None), help="Specify whether to use locally tools or globally installed tools.", hidden=False, nargs=1)
@click.option("--guess-ftype", type=(bool), default=(None), help="Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=1)
@click.option("--mmoreseqs-ftype", type=(str, str, str, str), default=(None, None, None, None), help="Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=4)
@click.option("--alpha", type=(bool), default=(None), help="Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.", hidden=False, nargs=1)
@click.option("--beta", type=(bool), default=(None), help="Set MMORE beta parameter. Determines X-drop below global max score accepted.", hidden=False, nargs=1)
@click.option("--gamma", type=(bool), default=(None), help="Set MMORE gamma parameter. Determines X-drop below global max score before ending search.", hidden=False, nargs=1)
@click.option("--hard-limit", type=(bool), default=(None), help="Set MMORE hard limit parameter. Determines lowest permissable score before ending search.", hidden=False, nargs=1)
@click.option("--run-prep", type=(bool), default=(None), help="Run file preparation for the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-bias", type=(bool), default=(None), help="Use null bias during the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-full", type=(bool), default=(None), help="Run full quadratic search during the MMORE stage of pipeline. Overrides alpha/beta/gamma parameters.", hidden=False, nargs=1)
@click.option("--run-domains", type=(bool), default=(None), help="Run over all domains found by MMSEQS for all target-query pairs during MMORE stage of pipeline. Alternatively, only uses the searches the top domain.", hidden=False, nargs=1)
@click.option("--run-vit-mmore", type=(bool), default=(None), help="Run full viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit", type=(bool), default=(None), help="Run viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-post", type=(bool), default=(None), help="Run posterior during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--range", type=(bool), default=(None), help="Specify [0] start and [1] stop range of targets to search against the query database.", hidden=False, nargs=1)
@click.option("--run-filter", type=(bool), default=(None), help="Run all filters during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit-filter", type=(bool), default=(None), help="Run viterbi filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-cld-filter", type=(bool), default=(None), help="Run cloud filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-fwd-filter", type=(bool), default=(None), help="Run forward filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--vit-filter", type=(bool), default=(None), help="Set viterbi filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--cld-filter", type=(bool), default=(None), help="Set cloud filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--fwd-filter", type=(bool), default=(None), help="Set forward filter threshold cutoff score.", hidden=False, nargs=1)
def my_func(*args, **kwargs):
    prep_search(args, kwargs)
    pass


@cli.command(name="mmseqs-search", help="Runs MMSEQS stage of MMORESEQS search on target and query databases. Takes as arguments: [0] search target filepath, and [1] search query filepath.  Target filepath types: MMDB.  Query filepath types: MMDB.")
@click.argument("target_filepath", type=str)
@click.argument("query_filepath", type=str)
@click.option("--dbg", type=(bool), default=(None), help="Outputs debugger info.", hidden=True, nargs=1)
@click.option("--dbg-viz", type=(bool), default=(None), help="Outputs visualization debugger info.", hidden=True, nargs=1)
@click.option("--debugout", type=(bool), default=(None), help="Specify path to output debugger info.", hidden=True, nargs=1)
@click.option("--verbose", type=(int), default=(None), help="Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info", hidden=False, nargs=1)
@click.option("--num-threads", type=(int), default=(None), help="Number of parallel threads.", hidden=False, nargs=1)
@click.option("--search-type", type=(str), default=(None), help="Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence", hidden=False, nargs=1)
@click.option("--use-pvals", type=(bool), default=(None), help="Use P-values (as opposed to default E-values) for reporting and scoring thresholds.", hidden=False, nargs=1)
@click.option("--eval", type=(bool), default=(None), help="Set E-value filter threshold cutoff score for reporting.", hidden=False, nargs=1)
@click.option("--program-mmseqs", type=(str), default=(None), help="Specify location of MMSEQS program.", hidden=False, nargs=1)
@click.option("--program-hmmer", type=(str), default=(None), help="Specify location of HMMER program.", hidden=False, nargs=1)
@click.option("--program-mmoreseqs", type=(str), default=(None), help="Specify location of MMORESEQS program.", hidden=False, nargs=1)
@click.option("--script-dir", type=(str), default=(None), help="Specify location of MMORESEQS scripts directory.", hidden=False, nargs=1)
@click.option("--index", type=(str, str), default=(None, None), help="Specify location of [0] target and [1] query index files.", hidden=False, nargs=2)
@click.option("--local-tools", type=(bool), default=(None), help="Specify whether to use locally tools or globally installed tools.", hidden=False, nargs=1)
@click.option("--guess-ftype", type=(bool), default=(None), help="Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=1)
@click.option("--mmoreseqs-ftype", type=(str, str, str, str), default=(None, None, None, None), help="Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=4)
@click.option("--mmoreseqs-main-ftype", type=(bool), default=(None), help="Specify filetype of [0] target and [1] query input.", hidden=False, nargs=1)
@click.option("--temp", type=(bool), default=(None), help="Specify temporary working folder location.", hidden=False, nargs=1)
@click.option("--mmseqs-m8", type=(bool), default=(None), help="Specify location of input MMSEQS .m8 results file.", hidden=False, nargs=1)
@click.option("--dbsizes", type=(bool), default=(None), help="Specify size of [0] target and [1] query databases.", hidden=False, nargs=1)
def my_func(*args, **kwargs):
    mmseqs_search(args, kwargs)
    pass


@cli.command(name="mmore-search", help="Runs MMORE stage of MMORESEQS search on target and query databases, using results from MMSEQS stage. Takes as arguments: [0] search target filepath, [1] search query filepath, and [2] temporary working directory path. Target filepath types: HMM. Query filepath types: FASTA.")
@click.argument("target_filepath", type=str)
@click.argument("query_filepath", type=str)
@click.option("--dbg", type=(bool), default=(None), help="Outputs debugger info.", hidden=True, nargs=1)
@click.option("--dbg-viz", type=(bool), default=(None), help="Outputs visualization debugger info.", hidden=True, nargs=1)
@click.option("--debugout", type=(bool), default=(None), help="Specify path to output debugger info.", hidden=True, nargs=1)
@click.option("--verbose", type=(int), default=(None), help="Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info", hidden=False, nargs=1)
@click.option("--num-threads", type=(int), default=(None), help="Number of parallel threads.", hidden=False, nargs=1)
@click.option("--search-type", type=(str), default=(None), help="Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence", hidden=False, nargs=1)
@click.option("--use-pvals", type=(bool), default=(None), help="Use P-values (as opposed to default E-values) for reporting and scoring thresholds.", hidden=False, nargs=1)
@click.option("--eval", type=(bool), default=(None), help="Set E-value filter threshold cutoff score for reporting.", hidden=False, nargs=1)
@click.option("--program-mmseqs", type=(str), default=(None), help="Specify location of MMSEQS program.", hidden=False, nargs=1)
@click.option("--program-hmmer", type=(str), default=(None), help="Specify location of HMMER program.", hidden=False, nargs=1)
@click.option("--program-mmoreseqs", type=(str), default=(None), help="Specify location of MMORESEQS program.", hidden=False, nargs=1)
@click.option("--script-dir", type=(str), default=(None), help="Specify location of MMORESEQS scripts directory.", hidden=False, nargs=1)
@click.option("--index", type=(str, str), default=(None, None), help="Specify location of [0] target and [1] query index files.", hidden=False, nargs=2)
@click.option("--local-tools", type=(bool), default=(None), help="Specify whether to use locally tools or globally installed tools.", hidden=False, nargs=1)
@click.option("--guess-ftype", type=(bool), default=(None), help="Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=1)
@click.option("--mmoreseqs-ftype", type=(str, str, str, str), default=(None, None, None, None), help="Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=4)
@click.option("--mmoreseqs-main-ftype", type=(bool), default=(None), help="Specify filetype of [0] target and [1] query input.", hidden=False, nargs=1)
@click.option("--temp", type=(bool), default=(None), help="Specify temporary working folder location.", hidden=False, nargs=1)
@click.option("--mmseqs-m8", type=(bool), default=(None), help="Specify location of input MMSEQS .m8 results file.", hidden=False, nargs=1)
@click.option("--dbsizes", type=(bool), default=(None), help="Specify size of [0] target and [1] query databases.", hidden=False, nargs=1)
@click.option("--alpha", type=(bool), default=(None), help="Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.", hidden=False, nargs=1)
@click.option("--beta", type=(bool), default=(None), help="Set MMORE beta parameter. Determines X-drop below global max score accepted.", hidden=False, nargs=1)
@click.option("--gamma", type=(bool), default=(None), help="Set MMORE gamma parameter. Determines X-drop below global max score before ending search.", hidden=False, nargs=1)
@click.option("--hard-limit", type=(bool), default=(None), help="Set MMORE hard limit parameter. Determines lowest permissable score before ending search.", hidden=False, nargs=1)
@click.option("--run-prep", type=(bool), default=(None), help="Run file preparation for the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-bias", type=(bool), default=(None), help="Use null bias during the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-full", type=(bool), default=(None), help="Run full quadratic search during the MMORE stage of pipeline. Overrides alpha/beta/gamma parameters.", hidden=False, nargs=1)
@click.option("--run-domains", type=(bool), default=(None), help="Run over all domains found by MMSEQS for all target-query pairs during MMORE stage of pipeline. Alternatively, only uses the searches the top domain.", hidden=False, nargs=1)
@click.option("--run-vit-mmore", type=(bool), default=(None), help="Run full viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit", type=(bool), default=(None), help="Run viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-post", type=(bool), default=(None), help="Run posterior during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--range", type=(bool), default=(None), help="Specify [0] start and [1] stop range of targets to search against the query database.", hidden=False, nargs=1)
@click.option("--run-filter", type=(bool), default=(None), help="Run all filters during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit-filter", type=(bool), default=(None), help="Run viterbi filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-cld-filter", type=(bool), default=(None), help="Run cloud filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-fwd-filter", type=(bool), default=(None), help="Run forward filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--vit-filter", type=(bool), default=(None), help="Set viterbi filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--cld-filter", type=(bool), default=(None), help="Set cloud filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--fwd-filter", type=(bool), default=(None), help="Set forward filter threshold cutoff score.", hidden=False, nargs=1)
def my_func(*args, **kwargs):
    mmore_search(args, kwargs)
    pass


@cli.command(name="easy-search", help="Prepares file and completes full search on target and query databases. Takes as arguments: [0] search target filepath, [1] search query filepath, and [2] temporary working directory path. Target filepath types: HMM. Query filepath types: HMM, FASTA.")
@click.argument("target_filepath", type=str)
@click.argument("query_filepath", type=str)
@click.argument("mmseqs_m8_filepath", type=str)
@click.option("--dbg", type=(bool), default=(None), help="Outputs debugger info.", hidden=True, nargs=1)
@click.option("--dbg-viz", type=(bool), default=(None), help="Outputs visualization debugger info.", hidden=True, nargs=1)
@click.option("--debugout", type=(bool), default=(None), help="Specify path to output debugger info.", hidden=True, nargs=1)
@click.option("--verbose", type=(int), default=(None), help="Amount of output: 0=quiet, 1=errors, 2=warnings, 3=info", hidden=False, nargs=1)
@click.option("--num-threads", type=(int), default=(None), help="Number of parallel threads.", hidden=False, nargs=1)
@click.option("--search-type", type=(str), default=(None), help="Set type of search: P2S=profile-to-sequence, S2S=sequence-to-sequence", hidden=False, nargs=1)
@click.option("--use-pvals", type=(bool), default=(None), help="Use P-values (as opposed to default E-values) for reporting and scoring thresholds.", hidden=False, nargs=1)
@click.option("--eval", type=(bool), default=(None), help="Set E-value filter threshold cutoff score for reporting.", hidden=False, nargs=1)
@click.option("--program-mmseqs", type=(str), default=(None), help="Specify location of MMSEQS program.", hidden=False, nargs=1)
@click.option("--program-hmmer", type=(str), default=(None), help="Specify location of HMMER program.", hidden=False, nargs=1)
@click.option("--program-mmoreseqs", type=(str), default=(None), help="Specify location of MMORESEQS program.", hidden=False, nargs=1)
@click.option("--script-dir", type=(str), default=(None), help="Specify location of MMORESEQS scripts directory.", hidden=False, nargs=1)
@click.option("--local-tools", type=(bool), default=(None), help="Specify whether to use locally tools or globally installed tools.", hidden=False, nargs=1)
@click.option("--guess-ftype", type=(bool), default=(None), help="Specify whether to infer filetype of target and query input. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=1)
@click.option("--mmoreseqs-ftype", type=(str, str, str, str), default=(None, None, None, None), help="Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query.  Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.", hidden=False, nargs=4)
@click.option("--alpha", type=(bool), default=(None), help="Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.", hidden=False, nargs=1)
@click.option("--beta", type=(bool), default=(None), help="Set MMORE beta parameter. Determines X-drop below global max score accepted.", hidden=False, nargs=1)
@click.option("--gamma", type=(bool), default=(None), help="Set MMORE gamma parameter. Determines X-drop below global max score before ending search.", hidden=False, nargs=1)
@click.option("--hard-limit", type=(bool), default=(None), help="Set MMORE hard limit parameter. Determines lowest permissable score before ending search.", hidden=False, nargs=1)
@click.option("--run-prep", type=(bool), default=(None), help="Run file preparation for the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-bias", type=(bool), default=(None), help="Use null bias during the MMORE stage of pipeline", hidden=False, nargs=1)
@click.option("--run-full", type=(bool), default=(None), help="Run full quadratic search during the MMORE stage of pipeline. Overrides alpha/beta/gamma parameters.", hidden=False, nargs=1)
@click.option("--run-domains", type=(bool), default=(None), help="Run over all domains found by MMSEQS for all target-query pairs during MMORE stage of pipeline. Alternatively, only uses the searches the top domain.", hidden=False, nargs=1)
@click.option("--run-vit-mmore", type=(bool), default=(None), help="Run full viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit", type=(bool), default=(None), help="Run viterbi during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-post", type=(bool), default=(None), help="Run posterior during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--range", type=(bool), default=(None), help="Specify [0] start and [1] stop range of targets to search against the query database.", hidden=False, nargs=1)
@click.option("--run-filter", type=(bool), default=(None), help="Run all filters during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-vit-filter", type=(bool), default=(None), help="Run viterbi filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-cld-filter", type=(bool), default=(None), help="Run cloud filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--run-fwd-filter", type=(bool), default=(None), help="Run forward filter during the MMORE stage of pipeline.", hidden=False, nargs=1)
@click.option("--vit-filter", type=(bool), default=(None), help="Set viterbi filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--cld-filter", type=(bool), default=(None), help="Set cloud filter threshold cutoff score.", hidden=False, nargs=1)
@click.option("--fwd-filter", type=(bool), default=(None), help="Set forward filter threshold cutoff score.", hidden=False, nargs=1)
def my_func(*args, **kwargs):
    easy_search(args, kwargs)
    pass


@cli.command(name="version", help="Get version of MMOREseqs.")
def my_func(*args, **kwargs):
    mmoreseqs_version(args, kwargs)
    pass


### COMMAND HELPER FUNCTIONS ###

def build_commandline_args(cmd_name, *args, **kwargs):
    arg_dict = args[0][1]
    opt_str = ""
    arg_str = ""
    for key in arg_dict.keys():
        if (key in primary_args[cmd_name]):
            arg_str += f"{arg_dict[key]} "
            continue
        if (type(arg_dict[key]) == tuple):
            value = arg_dict[key][0]
            if (value != None):
                opt_str += f"--{key} "
                for i in range(len(arg_dict[key])):
                    opt_str += f"--{arg_dict[key][i]} "
        else:
            value = arg_dict[key]
            if (arg_dict[key] != None):
                opt_str += f" --{key} {arg_dict[key]} "
    cmd_str = f"{cmd_name} {arg_str} {opt_str}"
    return cmd_str


def run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs):
    commandline_str = build_commandline_args(cmd_name, args, kwargs)
    print(f"RUN_COMMAND: {commandline_str}")
    # mmoreseqs_pylib.run(commandline_str)
    subprocess.run(f"mmoreseqs {commandline_str}")
    pass

### COMMAND FUNCTIONS ###


def prep(*args, **kwargs):
    cmd_name = "prep"
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


def prep_search(*args, **kwargs):
    cmd_name = "prep-search"
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


def easy_search(*args, **kwargs):
    cmd_name = "easy-search"
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


def mmseqs_search(*args, **kwargs):
    cmd_name = "mmseqs-search"
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


def mmore_search(*args, **kwargs):
    cmd_name = "mmore-search"
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


def mmoreseqs_version(*args, **kwargs):
    cmd_name = "version"
    # mmoreseqs_pylib.get_version()
    run_mmoreseqs_from_commandline_args(cmd_name, *args, **kwargs)
    pass


### MAIN ###

if __name__ == '__main__':
    cli()
