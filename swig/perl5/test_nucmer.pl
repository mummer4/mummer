use strict;
use warnings;
use Test::More;

require_ok('mummer');

my @bases = qw(A C G T);
sub seq {
  my ($n) = @_;
  my $s;
  for(my $i = 0; $i < $n; ++$i) {
    $s .= $bases[int(rand(4))];
  }
  return $s;
}

# Exact pair
{
  my $s1 = seq(100);
  my $s2 = substr($s1, 80) . seq(80);
  my $o = mummer::Options->new;
  $o->minmatch(10);
  $o->mincluster(10);
  my $a = mummer::align_sequences($s1, $s2, $o);

  ok(@$a == 1, "nb matches");
  ok($$a[0]{dirB} == 1, "direction of match");
  ok($$a[0]{sA} == 81, "sA");
  ok($$a[0]{eA} == 100, "eA");
  ok($$a[0]{sB} == 1, "sB");
  ok($$a[0]{eB} == 20, "eB");
  ok($$a[0]{Errors} == 0, "Errors");
  ok($$a[0]{SimErrors} == 0, "SimErrors");
  ok($$a[0]->identity == 1.0, "Identity");
  ok($$a[0]->similarity == 1.0, "Similarity");
}

done_testing;
