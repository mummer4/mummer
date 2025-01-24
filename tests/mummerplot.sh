mummer -mum -b -c -n $D/small_genome.fa $D/small_reads_0.fa > small_0.mums

mummerplot -t png -p with_spaces small_0.mums
sed 's/^ \+//' < small_0.mums > small_0_nospace.mums
mummerplot -t png -p no_space small_0_nospace.mums

diff with_spaces.fplot no_space.fplot
diff with_spaces.rplot no_space.rplot
