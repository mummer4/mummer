time nucmer --mum --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 81fe30ba229be8f0014d754c39cd3287

time nucmer --maxmatch --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 fa61620d01b700f476b6a19d3af28056

time nucmer --maxmatch --delta /dev/stdout -L 90 $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 ff9433627943d8ededbb70dcfa80b3dd

time nucmer --mum --large --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 81fe30ba229be8f0014d754c39cd3287

time nucmer --maxmatch --large --delta /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 fa61620d01b700f476b6a19d3af28056

time nucmer --maxmatch --large --delta /dev/stdout -L 90 $D/seed_reads_1.fa $D/seed_reads_0.fa | \
    ufasta sort -H | test_md5 ff9433627943d8ededbb70dcfa80b3dd
