# MMOREseqs

Pruned Forward-Backward implementation of profile HMM alignment.

## Installation / Getting Started

### Build from Source

We use a "git flow" workflow. We have one active branch:
* **master** will be the stable MMOREseqs release branch.  

To clone your own copy of the MMOREseqs repository for the first time:

```bash
   $ git clone https://github.com/TravisWheelerLab/mmoreseqs mmoreseqs
   $ cd mmoreseqs
   $ make
```

Make will build an executable binary and a python module in the `mmoreseqs/build\bin` directory.  

For more information about gitflow, see the
[mmoreseqs wiki](https://github.com/TravisWheelerLab/mmoreseqs/wiki)

The executable, called `mmoreseqs` will end up in the `build/` subdirectory.
Dependencies, specifically [Easel](https://github.com/EddyRivasLab/easel) will
be fetched automatically at the correct versions.

### Example 

For a quick example database to learn the workflow and check the correct behavior.  This test is contained in `mmoreseqs/example`.  In this file, we have a collection of shell scripts and small target and query databases  in `example/db`.  There are a few workflow examples shown: 
- (1) The first, `easysearch-example.sh`.  Using easy-search allows the user to search the databases with one single command, just specify the target and query, as well as the temporary working directory to use.
- (2a) The second, `prepsearch-example-1.sh` and `prepsearch-example-2.sh`.  Using this workflow splits the preparation step from the search step.  This can bypass alot of unneccessary overhead, especially when performing multiple searches against a particular database under different parameters.
- (2b) The third breaks down `prepsearch-example-2.sh` further breaks the search down into `prepsearch-example-2-mmseqs.sh`, the first MMseqs2 search phase and `prepsearch-example-2-mmore.sh`, the second MMORE search phase.  All of these methods use the same parameters and will generate the same results.

## Usage

### Workflow Pipelines

The MMOREseqs workflow takes in three primary inputs:

- `QUERY`
  - A multiple sequence alignment (MSA) database file.  
- `TARGET`
  - A sequence (FASTA) database file that the query will search against.
- `TEMP_DIR`
  - A working directory for storing intermediate data during the search workflows. Can

There are a few ways to execute the MMOREseqs Pipeline.  

(1) Simplest method: `mmoreseqs easy-search`
Easy-search handles all the input files through the entire search pipeline in a single command.  

```
mmoreseqs easy-search <i:query_msa> <i:target_fasta> <i:temp_dir>
```

- Arguments:
  - `<query_msa>`
    - Query multiple sequence alignment (MSA) database.  
  - `<target_fasta>`
    - Target sequence (FASTA) database to be searched against.  
  - `<temp_dir>`
    - A temporary working directory. Pipeline will create this directory so long as its parent directory exists.  
  

(2) Next simplest method: `mmoreseqs prep-search`
Prep-search runs the search workflow in two stages: `prep`, which prepares the databases in a format in a way in which it can be easily consumed, and; `prep-search` which performs the search given the expected prep directory.

```
mmoreseqs prep <i:query_msa> <i:target_fasta> <i:db_dir>
mmoreseqs prep-search <i:db_dir>
```

- Arguments:
  - `<query_msa>`
    - Query multiple sequence alignment (MSA) database.  
  - `<target_fasta>`
    - Target sequence (FASTA) database to be searched against.  
  - `<temp_dir>`
    - A temporary working directory. Pipeline will create this directory so long as its parent directory exists.  


(3) Highest control method: `mmoreseqs search`
Search runs the main MMseqs and MMORE steps of the MMOREseqs workflow.  

```
mmoreseqs search <i:query_mmore_hmm> <i:target_mmore_fasta> \
   <i:query_mmseqs_smmdb> <i:target_mmseqs_pmmdb> <i:target_mmseqs_smmdb>
```

- Arguments:
  - `<query_mmore_hmm>`
    - Query Hidden Markov Model (HMM) database file. HMMER format.  
  - `<target_mmore_fasta>`
    - Target sequence (FASTA) database file.  
  - `<query_mmseqs_smmdb>`
    - Query sequence (SMMDB) database file. MMSEQS format.  
  - `<query_mmseqs_pmmdb>`
    - Query profile (PMMDB) database file. MMSEQS format.  
  - `<target_mmseqs_smmdb>`
    - Target sequence (SMMDB) database file. MMSEQS format.  


(4) Subroutines: `mmoreseqs mmore-search` `mmoreseqs mmseqs-search`  
Independently searches the MMseqs or the MMORE steps of MMOREseqs workflow. Because of the input/output dependencies, remember that the MMseqs stage must always be performed before the MMORE stage..

```
mmoreseqs mmseqs-search <i:query_mmseqs_pmmdb> <i:query_mmseqs_smmdb> <i:target_mmseqs_smmdb>
```

- Arguments:
  - `<query_mmseqs_smmdb>`
    - Query profile (PMMDB) database file. MMSEQS format.  
  - `<query_mmseqs_pmmdb>`
    - Query sequence (SMMDB) database file. MMSEQS format.  
  - `<target_mmseqs_smmdb>`
    - Target sequence (SMMDB) database file. MMSEQS format. 

```
mmoreseqs mmore-search <i:query_mmore_hmm> <i:target_mmore_fasta> <i:results_mmseqs_m8>
```

- Arguments:
  - `<query_mmore_hmm>`
    - Query profile (HMM) database file. HMMER format.
  - `<target_mmore_fasta>`
    - Target sequence (FASTA) database file.
  - `<results_mmseqs_m8>`
    - Results file (.m8) outputted from MMseqs stage of pipeline.

### Workflow Options

- General Options:
  - `--verbose INT`
    - The amount of output: 0=quiet, 1=errors, 2=warnings, 3=info
  - `--num-threads INT`
    - The number of parallel threads to run.  Does not currently work with all workflows. 
  - `--eval DOUBLE`
    - Set E-value filter threshold cutoff score for reporting.
  - `--use-pvals BOOL`
    - Whether to use P-value (as opposed to using the default E-value) for assessing reporting and filtering thresholds. Uses value stored in `--eval` for reporting threshold.

- Pipeline Options (Used with `mmoreseqs search`.):
   `--run-prep BOOL`
    - Run file preparation stage of pipeline.
  - `--run-mmseqs BOOL`
    - Run MMseqs stage of the pipeline.
  - `--run-mmseqs-pref BOOL`
    - Run MMseqs prefilter stage of the pipeline.
  - `--run-mmseqs-align BOOL`
    - Run MMseqs alignment stage of the pipeline.
  - `--run-convert BOOL`
    - Run conversion of MMseqs output to MMORE input stage of the pipeline.
  - `--run-mmore BOOL`
  - - Run MMORE stage of the pipeline.
  - `--run-vit BOOL`
    - Run viterbi during the MMORE stage of pipeline. Otherwise bypasses viterbi computation.
  - `--run-post BOOL`
    - Run posterior during the MMORE stage of pipeline. Otherwise bypasses posterior computation.

- Input File Options:
  - `--tmp/--prep`
    - Set temporary preparation working directory for pipeline interrim operations.
  - `--mmseqs-m8 TEXT`
    - Set output .m8 file from MMseqs search as input for MMORE search.
  - `--index TEXT TEXT`  
    - Set index file of [0] query HMM and [1] target FASTA file for faster access.

- MMseqs Options (More information can be found about these options in the MMseqs2 User Guide):
  - `--mmseqs-kmer INT=7`
    - Sets the kmer length during the prefilter step of MMseqs stage.
  - `--mmseqs-kscore INT=80`
    - Sets filter score for kmers during the prefilter step of MMseqs stage.
- - `--mmseqs-hits-per-search INT=1000`
    - Sets total number of hits per target allowed per search to be reported by prefilter step of MMseqs stage.
  - `--mmseqs-ungapped-vit INT=15`
    - Sets filter score for ungapped viterbi step of MMseqs stage.
  - `--mmseqs-eval FLOAT=1e-2`
    - Sets filter E-value of MMseqs for (gapped) viterbi step of MMseqs stage.
  - `--mmseqs-pval FLOAT`
    - Sets filter P-value of MMseqs for (gapped) viterbi step of MMseqs stage. Overrides E-value score.
  - `--mmseqs-altalis INT=0`
    - Sets number of alternate alignments that can be found per 
  - `--mmseqs-split INT`
    - When building MMseqs databases determines the number of independent files to store databases across. Larger number keeps file size smaller.

- MMORE Options:
  - `--vit-filter FLOAT=1e-2` 
    - Set viterbi filter threshold cutoff score.
  - `--cld-filter FLOAT=1e-3` 
    - Set cloud filter threshold cutoff score.
  - `--fwd-filter FLOAT=1e-4`
    - Set bounded forward filter threshold cutoff score.
  - `--alpha FLOAT=12.0`
    - Set MMORE alpha parameter. Determines X-drop below current anti-diagonal max score accepted.
  - `--beta FLOAT=16.0`
    - Set MMORE beta parameter. Determines X-drop below global max score accepted before termination of search.
  - `--gamma INT=5`
    - Set MMORE gamma parameter. Determines number of anti-diagonal sweeps performed before beginning pruning.
  - `--hard-limit FLOAT=(-12.0)`
    - Set MMORE hard limit parameter. Determines lowest permissable score before termination of search.
  - `--range INT INT` 
    - Specify [0] start and [1] stop range of .m8 MMseqs results to search.

- Output File Options (These specify types of output.):
  - `--stdout TEXT`
    - Redirects standard output to specified file.
  - `--stderr TEXT`
    - Redirects standard error to specified file.
  - `--m8out TEXT`
    - Outputs the results of MMOREseqs search to MMseqs-style .m8 format tsv file.
  - `--domtblout TEXT`
    - Outputs the results of MMOREseqs search to HMMER-style domain table format tsv file.
  - `--mmseqs-m8out TEXT`
    - Outputs the results of MMseqs stage to .m8 format file.

- Various Uncommon Options (Users should rarely need to invoke these except in odd cases.):
  - `--run-domains BOOL`
    - Run over all domains found by MMSEQS for all target/query pairs during MMORE stage of pipeline. Alternatively, this only searches the top domain.
  - `--run-bias BOOL`
    - Use null bias during the MMORE stage of pipeline.
  - `--run-filter BOOL`
    - Run all filters during the MMORE stage of pipeline.
  - `--run-vit-filter BOOL` 
    - Run viterbi filter during the MMORE stage of pipeline. Still computes viterbi score.
  - `--run-cld-filter BOOL` 
    - Run cloud filter during the MMORE stage of pipeline. Still computes cloud score.
  - `--run-fwd-filter BOOL`
    - Run forward filter during the MMORE stage of pipeline. Still computes forward score.
  - `--run-vit-mmore BOOL`
    - Run full viterbi during the MMORE stage of pipeline. Bypasses viterbi computation.
  - `--run-full BOOL`
    - Run full quadratic search during the MMORE stage of pipeline.Overrides alpha/beta/gamma parameters.
  
- Forwarding Options (General users should never need to invoke these options. They store data which are passed to scripts for intermediary workflow operations.):
  - `--search-type TEXT`
    - Set type of search inputs. P2S=profile-to-sequence search, S2S=sequence-to-sequence search.
  - `--guess-ftype BOOL`
    - Specify whether to infer filetype of target and query input by reading file extensions. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.
  - `--mmoreseqs-ftype <TEXT TEXT TEXT TEXT>`
    - Specify filetype of [0] target, [1] query, [2] MMSEQS target, [3] MMSEQS query. Options: HMM, FASTA, MSA, MM_MSA, HHM, MMDB, MMDB_S, MMDB_P.
  - `--program-mmseqs TEXT`
    - Specify location of MMSEQS program.
  - `--program-hmmer TEXT`
    - Specify location of HMMER program.
  - `--program-mmoreseqs TEXT`
    - Specify location of MMORESEQS program.
  - `--script-dir TEXT`
    - Specify location of MMORESEQS scripts directory.
  - `--local-tools BOOL`
    - Specify whether to use locally tools or globally installed tools.
  - `--dbsizes INT INT`
    - Specify query and target database sizes. Used when database computed E-values on a split database.

- Developer Options (General users whould never need to invoke these options.  Used by developers for testing and debugging. Many require that project is built with `$ make build-debug`.):
  - `--dbg BOOL`
    - Whether to output debugger information.
  - `--dbg-viz BOOL`
    - Whether to output debugger visualization information.
  - `--allout TEXT`
    - Outputs a exhaustive tsv file.
  - `--myout TEXT`
    - Outputs a custom formatted tsv file.
  - `--mytimeout TEXT`
    - Outputs a runtime tsv file.
  - `--mythreshout TEXT`
    - Outputs a threshold score tsv file.
  - `--myhmmerout TEXT`
    - Outputs a HMMER-style file.
  - `--debugout TEXT`
    - Output debugger information to file.

## Development

Python dependencies can be installed in virtual environment with `% make create-env`.

Build-Time Dependencies:

- easel (https://github.com/EddyRivasLab/easel)
- pybind11 (v2.8.1)

Run-Time Dependencies:

- HMMER (v3.3.2) (http://hmmer.org/download.html)
- MMseqs2 (https://github.com/soedinglab/MMseqs2.git)
   - commit hash: fa4cd2a7ab15a2636b5c1615859a2491d6888300  
- Python (v2.7)
