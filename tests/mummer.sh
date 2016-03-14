mummer -mum $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta hsort -H | ufasta dsort | test_md5 4e8182c9f745abf59158f69a05b942f3
mummer -maxmatch $D/seed_reads_1.fa $D/seed_reads_0.fa | ufasta hsort -H | ufasta dsort | test_md5 a459f93742d1c36819e53e7a4c128bf7
