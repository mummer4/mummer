# SWIG binding installation

## Dependencies

If compiling from a [release source tarball](releases), the additional
dependencies are the source headers and library of the target scripting
language (Python, Ruby or Perl).

If compiling from the gihub development tree, the swig wrapper files
are not included. You need [SWIG](http://www.swig.org/) version >=
3.0 to generate the wrappers.

## Configuration

In the following, we use Python as the default target scripting
languages. The instructions are similar for other language by
replacing `python` with `perl` or `ruby`.

Compilation from a release tarball:

```Shell
./configure --enable-python-binding
```

Compilation from the github development tree:

```Shell
./configure --enable-swig --enable-python-binding
```

## Installation directory

There are 3 possibilities on where the binding will be installed:

1. If no `--prefix` option is passed to `./configure`, the library
will be installed in the normal system directories, as defined at
compilation time of Python (on my system that would be
`/usr/lib/python2.7/dist-packages`).

2. If `--prefix` is used, then the binding is installed in
`$(prefix)/lib/python`.

3. A directory for installation can be passed to `--enable-python-binding`, like so:

```Shell
./configure --enable-python-binding=$HOME/lib/python
```

## Using the bindings

Unless the bindings are installed in the default system directories
(option 1 above), one must set environment variables for the scripting
language to find the binding:

* `PYTHONPATH` for Python
* `PERL5LIB` for Perl
* `RUBYLIB` for Ruby

If the installation is successful and the environment is properly set,
then the following should print OK:

```Shell
python -c 'import mummer; print("OK")'
ruby -rmummer -e 'puts("OK")'
perl -Mmummer -e 'print("OK\n")'
```
