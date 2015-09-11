set -e
set -x

../nucmer --save small_sa --delta small_1.delta small_reads_1.fa small_reads_0.fa
../nucmer --load small_sa --delta small_2.delta small_reads_1.fa small_reads_0.fa
diff <(ufasta sort -H small_1.delta) <(ufasta sort -H small_2.delta)
