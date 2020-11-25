# fb-pruner

Pruned Forward-Backward implementation of profile HMM alignment

## Development

### Setup

This project uses Autotools. To create the script that will configure a Makefile: 
`autoreconf -i`

If building from source, must pull in git submodule:
`git submodule update --init --recursive`
or 
`git submodule update --recursive --remote`

### Build

To build from source:

The usual incantation will work:

```
autoreconf -i
./configure
make
```

This will produce a binary called `mmore_seqs`.

### Usage

```
mmore_seqs
```