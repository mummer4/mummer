#! /usr/bin/env perl

use strict;
use warnings;
use mummer;

# Read 2 files that are raw sequences: no fasta header or spaces or
# anything else but sequence. It aligns the two sequences with some
# default parameters and then output the results in a format
# comparable to the delta format.

sub read_file {
  local $/ = undef;
  my ($file) = @_;
  open(FILE, "<", $file) or die "Can't open file '$file': $!";
  binmode FILE;
  my $s = <FILE>;
  close(FILE);
  return $s;
}

my $ref = read_file($ARGV[0]);
my $qry = read_file($ARGV[1]);

my $o = mummer::Options->new;
$o->minmatch(10);
$o->mincluster(10);
my $aligns = mummer::align_sequences($ref, $qry, $o);
print($aligns, "\n");

for my $a (@$aligns) {
  print("$$a{sA} $$a{eA} $$a{sB} $$a{eB} $$a{Errors} $$a{SimErrors} $$a{NonAlphas}\n");
  # Need to print the deltas
  print($$a{delta}, "\n");
  while(my ($k, $v) = each %{$$a{delta}}) {
    print("$k $v\n");
  }
  print("0\n");
}
