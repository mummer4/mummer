essamummer -mum $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta dsort /dev/stdin | test_md5 66805ce38f31278c960ea28787494f02
essamummer -maxmatch $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta dsort /dev/stdin | test_md5 5ae11b38619c97ce2537f97ce90fb162
