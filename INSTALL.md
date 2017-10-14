## MUMmer4 INSTALLATION README

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

## Compilation & Installation

To compile and install:

```Shell
./configure --prefix=/path/to/installation
make
make install
```

If `--prefix` is omitted, the software is installed in
`/usr/local`. One may need `sudo make install` if installing in a
system location.

If compiling from the git tree, do `autoreconf -fi` first.

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



# UN-INSTALLATION

In the "MUMmer4.x" base directory type:

```Shell
make uninstall
```

## CONTACT INFORMATION

Please address questions and bug reports via the [github issue
tracker](../../issues).
