nucmer --delta batch.delta --maxmatch --batch 5000 $D/small_reads_1.fa $D/small_reads_0.fa
nucmer --delta no_batch.delta --maxmatch $D/small_reads_1.fa $D/small_reads_0.fa
diff <(ufasta sort -H no_batch.delta) <(ufasta sort -H batch.delta)
