nucmer --save ${N}_sa --delta ${N}_1.delta $D/small_reads_1.fa $D/small_reads_0.fa
nucmer --load ${N}_sa --delta ${N}_2.delta $D/small_reads_1.fa $D/small_reads_0.fa
diff <(ufasta sort -H ${N}_1.delta) <(ufasta sort -H ${N}_2.delta) > $N.diff
