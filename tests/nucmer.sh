# nucmer --mum --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta sort /dev/stdin | test_md5 43c892f87150a8503ac48ce7569369f7
nucmer --maxmatch --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta sort /dev/stdin | test_md5 7295d9568b9f14e3eb5ac0283937446f
