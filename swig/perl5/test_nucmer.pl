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

  ok(@$a >= 1, "nb matches");
  my @b = grep { $$_{sA} == 81 && $$_{eA} == 100 } @$a;
  ok(@b == 1, "one perfect match");
  ok($b[0]{dirB} == 1, "direction of match");
  ok($b[0]{sA} == 81, "sA");
  ok($b[0]{eA} == 100, "eA");
  ok($b[0]{sB} == 1, "sB");
  ok($b[0]{eB} == 20, "eB");
  ok($b[0]{Errors} == 0, "Errors");
  ok($b[0]{SimErrors} == 0, "SimErrors");
  ok($b[0]->identity == 1.0, "Identity");
  ok($b[0]->similarity == 1.0, "Similarity");
}

# Attempt changing number of threads
my $nt = mummer::get_num_threads();
ok($nt >= 1, "Number threads");
mummer::set_num_threads($nt + 1);
my $nnt = mummer::get_num_threads();
ok($nnt == 1 || $nnt == $nt + 1, "Increased threads");
print("$nt $nnt\n");

done_testing;
