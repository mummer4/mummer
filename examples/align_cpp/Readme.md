# Description

Similar to the scripts `align.*`, this program is an example to show
how to align two sequences. This program is very primitive and
meant only as an example:

- there is little error handling
- it takes two files containing sequences (no fasta format, only the
sequence)
- the parameters to the aligner are hard coded

# Compilation

Make sure that you installed mummer library properly. The following
command should print "OK": `pkg-config --exists mummer && echo OK`. If
not, check the followings:

- that you ran `make install`
- if you used `--prefix /path` with `./configure`, you may have to
append `/path/lib/pkgconfig` to the environment variable `PKG_CONFIG_PATH`.

Now `make` should create the executable `align`.
