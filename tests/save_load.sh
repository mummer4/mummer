nucmer --save ${N}_sa0 $D/small_reads_1.fa
nucmer --save ${N}_sa1 --delta ${N}_1.delta $D/small_reads_1.fa $D/small_reads_0.fa
nucmer --load ${N}_sa1 --delta ${N}_2.delta $D/small_reads_1.fa $D/small_reads_0.fa
nucmer --load ${N}_sa0 --delta ${N}_3.delta $D/small_reads_1.fa $D/small_reads_0.fa
diff <(ufasta sort -H ${N}_1.delta) <(ufasta sort -H ${N}_2.delta) > ${N}_2.diff
diff <(ufasta sort -H ${N}_1.delta) <(ufasta sort -H ${N}_3.delta) > ${N}_3.diff
