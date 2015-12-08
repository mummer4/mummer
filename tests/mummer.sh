mummer -mum $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta hsort -H | ufasta dsort | test_md5 dde330968ea5e0341a786f318ba5d6c9
mummer -maxmatch $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta hsort -H | ufasta dsort | test_md5 d75f4b0ca0a30cabc603ae8d30acf0f2
