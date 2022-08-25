# MMORESEQS

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

## Example 

For a quick example database to learn the workflow and check the correct behavior.  This test is contained in `mmoreseqs/example`.  In this file, we have a collection of shell scripts and small target and query databases  in `example/db`.  There are a few workflow examples shown: 
- (1) The first, `easysearch-example.sh`.  Using easy-search allows the user to search the databases with one single command, just specify the target and query, as well as the temporary working directory to use.
- (2a) The second, `prepsearch-example-1.sh` and `prepsearch-example-2.sh`.  Using this workflow splits the preparation step from the search step.  This can bypass alot of unneccessary overhead, especially when performing multiple searches against a particular database under different parameters.
- (2b) The third breaks down `prepsearch-example-2.sh` further breaks the search down into `prepsearch-example-2-mmseqs.sh`, the first MMseqs2 search phase and `prepsearch-example-2-mmore.sh`, the second MMORE search phase.  All of these methods use the same parameters and will generate the same results.

## Usage

MMOREseqs Workflow:  
The MMOREseqs takes in three primary inputs:

- QUERY: a multiple sequence alignment (MSA) database file.  
- TARGET: a sequence (FASTA) database file that the query will search against.

There are a few ways to execute the MMOREseqs Pipeline.  

(1) Simplest method: `mmoreseqs easy-search`
Easy-search handles all the input files through the entire search pipeline in a single command.  

Takes three arguments:

- query_msa: Query multiple sequence alignment (MSA) database.  
- target_fasta: Target sequence (FASTA) database to be searched against.  
- temp_dir: A temporary working directory. Pipeline will create this directory so long as its parent directory exists.  

```
mmoreseqs easy-search <i:query_msa> <i:target_fasta> <i:temp_dir>
```

(2) Next simplest method: `mmoreseqs prep-search`
Prep-search runs the search workflow in two stages: `prep`, which prepares the databases in a format in a way in which it can be easily consumed, and; `prep-search` which performs the search given the expected prep directory.

Arguments:

- query_msa: Query multiple sequence alignment (MSA) database.  
- target_fasta: Target sequence (FASTA) database to be searched against.  
- temp_dir: A temporary working directory. Pipeline will create this directory so long as its parent directory exists.  

```
mmoreseqs prep <i:query_msa> <i:target_fasta> <i:db_dir>
mmoreseqs prep-search <i:db_dir>
```

(3) Highest control method: `mmoreseqs search`
Search runs the main MMSEQS and MMORE steps of the MMORESEQS workflow.  

Arguments:

- query_mmore_hmm: Query Hidden Markov Model (HMM) database file. HMMER format.  
- target_mmore_fasta: Target sequence (FASTA) database file.  
- query_mmseqs_smmdb: Query sequence (SMMDB) database file. MMSEQS format.  
- query_mmseqs_pmmdb: Query profile (PMMDB) database file. MMSEQS format.  
- target_mmseqs_smmdb: Target sequence (SMMDB) database file. MMSEQS format.  

```
mmoreseqs search <i:query_mmore_hmm> <i:target_mmore_fasta> \
   <i:query_mmseqs_smmdb> <i:target_mmseqs_pmmdb> <i:target_mmseqs_smmdb>
```

(4) Subroutines: `mmoreseqs mmore-search` `mmoreseqs mmseqs-search`  

Arguments:

- query_mmore_hmm: Query Hidden Markov Model (HMM) database file. HMMER format.  
- target_mmore_fasta: Target sequence (FASTA) database file.  
- query_mmseqs_smmdb: Query sequence (SMMDB) database file. MMSEQS format.  
- query_mmseqs_pmmdb: Query profile (PMMDB) database file. MMSEQS format.  
- target_mmseqs_smmdb: Target sequence (SMMDB) database file. MMSEQS format.  

```
mmoreseqs mmseqs-search <i:query_mmseqs_pmmdb> <i:query_mmseqs_smmdb> <i:target_mmseqs_smmdb> <o:results_mmseqs_m8>
mmoreseqs mmore-search <i:query_mmore_hmm> <i:target_mmore_fasta> <i:results_mmseqs_m8>
```

## Development

Build-Time Dependencies:

- easel (https://github.com/EddyRivasLab/easel)

Run-Time Dependencies:

- HMMER (v3.3.2) (http://hmmer.org/download.html)
- MMseqs2 (https://github.com/soedinglab/MMseqs2.git)
   - commit hash: fa4cd2a7ab15a2636b5c1615859a2491d6888300  
- Python (v2.7)
