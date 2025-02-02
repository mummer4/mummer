# MUMmer4 Compilation README

## Dependencies

If compiling from a [release source tarball](../../releases) (recommended), you need a recent version of the GCC compiler (see below, only GCC is supported) and other essential tools for compilation (GNU make, ar, etc).
Additional requirements are needed to compile the SWIG script bindings.
See the [SWIG installation guide](swig/INSTALL.md).

If compiling from the github development tree, additionally you need autotools (autoconf, automake and libtools), [yaggo](https://github.com/gmarcais/yaggo/releases).

### On Ubuntu

From the tarball:
```Shell
sudo apt install build-essential
# For the bindings to scripting, additionally install
sudo apt install swig python3-dev ruby-dev libperl-dev
```

From the git tree:
```Shell
sudo apt instaoo build-essential git yaggo autoconf automake libtool gettext
# For the bindings to scripting, additionally install
sudo apt install swig python3-dev ruby-dev libperl-dev
```

### On Mac OS

MUMmer must be compiled with GCC, not Clang (and not the Apple provided `gcc` which is really `clang`).
Install with Brew:

```Shell
brew install autoconf automake libtool md5sha1sum
gem install yaggo
```

### On FreeBSD

MUMmer must be compiled with GCC, not Clang.
Install with Brew:

```Shell
brew install autoconf automake libtool md5sha1sum bash
gem install yaggo
```


## Compilation & Installation

If compiling from the release tarball (recommended), then the first command `autoreconf -fi` is not necessary.

### On Ubuntu

```Shell
autoreconf -fi # Optional, on if compiling from git tree
./configure --prefix=/path/to/installation
make
make check # Optional
make install
```

If `--prefix` is omitted, the software is installed in `/usr/local`.
One may need `sudo make install` if installing in a system location.

### On MacOS

Compile with `gcc-14`.

```Shell
autoreconf -fi # Optional, on if compiling from git tree
./configure --prefix=/path/to/installation CC=gcc-14 CXX=g++-14
make
make check # Optional
make install
```

### On FreeBSD

Compile with `gcc14`.

```Shell
autoreconf -fi # Optional, on if compiling from git tree
./configure --prefix=/path/to/installation MAKE=gmake CC=gcc14 CXX=g++14 LDFLAGS=-Wl,-rpath=/usr/local/lib/gcc14
gmake
gmake check # Optional
gmake install
```

## SOFTWARE REQUIREMENTS

The MUMmer4.x package requires the following to run successfully. In
the absence of one or more of these utilities, certain MUMmer programs
may fail to run correctly. In parenthesis the minimum version. These
utilities must be accessible via the system path:

* perl5 (5.6.0)
* sh
* sed
* awk

### OPTIONAL UTILITIES

To use the visualization tools included with MUMmer, it may be
necessary to install the following utilities:

* fig2dev (3.2.3)
* gnuplot (4.0)
* xfig    (3.2)



## UN-INSTALLATION

In the "MUMmer4.x" base directory type:

```Shell
make uninstall
```

## Container installation

### Docker container

To build a Docker container image containing mummer and all dependencies, [install Docker](https://docs.docker.com/get-docker/) for Windows, Mac, or Linux; clone the mummer git repo; and issue the following command in the top-level directory of the mummer git working tree:

```Shell
docker build -t mummer -f Dockerfiles/Dockerfile-debian .
```

To execute individual mummer commands (e.g., `nucmer`) within a container:

```Shell
docker run --rm -v $PWD:/mnt -w /mnt mummer nucmer -p <prefix> ref.fa  qry.fa
```

To execute an interactive shell within a container (from which mummer commands can be executed):

```Shell
docker run -it --rm -v $PWD:/mnt -w /mnt mummer
```

### Apptainer container

To build an Apptainer container, from the git tree run:

```Shell
apptainer build mummer.sif Apptainerfiles/mummer-alpine.def # or Apptainerfiles/mummer-debian.def
```

To build another branch than HEAD (for example `develop`), run:

```Shell
apptainer build --build-arg treeish=develop mummer.sif Apptainerfiles/mummer-alpine.def
```

To execute individual MUMmer commands (e.g., `nucmer`):

```Shell
apptainer run mummer.sif nucmer -p <prefix> ref.fa qry.fa
```

## CONTACT INFORMATION

Please address questions and bug reports via the [github issue
tracker](../../issues).
