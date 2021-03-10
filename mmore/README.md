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

New Install Pipeline:

```
cd {project_directory}/mmore/
make -f my_make/Makefile libs
make -f my_make/Makefile tools 
# (tools only necessary if hmmer and mmseqs are not installed on your system)
make -f my_make/Makefile 
make -f my_make/Makefile install INSTALL_DIR={install_directory} 
# (install_directory defaults to /bin/ in current directory)
```

This will produce a binary called `mmore`.

### Usage

```
mmore 
```