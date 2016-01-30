nucmer --mum --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 09a81068d91cb030a350ff5c8a999925

nucmer --maxmatch --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 d782ac6370f1d285828292ecbe8d7705

nucmer --maxmatch --delta /dev/stdout -L 90 $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 7f0dfdaf740f1a1889a28c0735d0ec46
