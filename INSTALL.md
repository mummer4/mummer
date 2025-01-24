# MUMmer4 INSTALLATION README

## Dependencies

If compiling from a [release source tarball](../../releases) you need a
recent version of the GCC compiler (g++ version >= 4.7) and other
essential tools for compilation (GNU make, ar, etc. Install
`build-essentials` on Debian or Ubuntu derivatives).  Additional
requirements are needed to compile the SWIG script bindings. See the
[SWIG installation guide](swig/INSTALL.md).

If compiling from the github development tree, additionally you need autotools (autoconf, automake and libtools),
[yaggo](https://github.com/gmarcais/yaggo/releases).
You should compile from a [release source tarball](../../releases), unless you plan on modifying the code of MUMmer.

On Ubuntu:

```Shell
sudo apt install git build-essential yaggo autoconf automake libtool gettext
# For the bindings to scripting, additionally install
sudo apt install swig python3-dev ruby-dev libperl-dev
```

## Compilation & Installation

To compile and install from a [release source tarball](../../releases):

```Shell
./configure --prefix=/path/to/installation
make
make install
```

If `--prefix` is omitted, the software is installed in
`/usr/local`. One may need `sudo make install` if installing in a
system location.

To compile from the github tree, `autoreconf` must additionally be run:
```Shell
autoreconf -fi
./configure --prefix=/path/to/installation
make
make install
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
docker build -t mummer .
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
apptainer build mummer.sif mummer.def
```

To build another branch than HEAD (for example `develop`), run:

```Shell
apptainer build --build-arg treeish=develop mummer.sif mummer.def
```

To execute individual MUMmer commands (e.g., `nucmer`):

```Shell
apptainer run mummer.sif nucmer -p <prefix> ref.fa qry.fa
```

## CONTACT INFORMATION

Please address questions and bug reports via the [github issue
tracker](../../issues).
