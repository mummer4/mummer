# Description

The scripts `align.pl`, `align.py` and `align.rb` are example scripts
to align two sequences and display the results in a format reminiscent
of the `.delta` format. These scripts are very primitive and meant
only as examples:

- no error handling
- they take two files containing sequences (no fasta format, only the
sequence)
- the parameters to the aligner are hard coded

# Known Limitations

1. The Python script works with Python2 but seems to crash on Python3
2. The Perl binding doe not return the position of the indels in the
alignment (delta array), only the general coordinates (SWIG problem?)
