nucmer --maxmatch --sam-long sam_test_long1.sam --save sam_test_index \
    $D/seed_reads_1.fa $D/seed_reads_0.fa
nucmer --maxmatch --sam-long sam_test_long2.sam --load sam_test_index \
    $D/seed_reads_1.fa $D/seed_reads_0.fa

pgline="$(printf "@PG\tID:nucmer\tPN:nucmer\tVN:%s\tCL:" "$VERSION")"
pglen=$(echo -n "$pgline" | wc -c)
for i in 1 2; do
    f=sam_test_long${i}.sam

    # Test cigar strings
    check_cigar $f $D/seed_reads_1.fa $D/seed_reads_0.fa

    # Test header
    [ "$(head -n 1 $f)" = "$(printf "@HD\tVN:1.4\tSO:unsorted")" ]
    [ "$(grep '^@PG' $f | head -c $pglen)" = "$pgline" ]

    # Test sequence headers
    diff -q <(grep '^@SQ' $f) <(ufasta sizes -H $D/seed_reads_1.fa | sed 's/\([0-9]\+\) \([0-9]\+\)/@SQ\tSN:\1\tLN:\2/')

    # Test that samtools can parse our output
    if [ -n "$SAMTOOLS" ]; then
        "$SAMTOOLS" view -b $f > /dev/null
    fi
done

# Test overall
nucmer --maxmatch --sam-short /dev/stdout $D/seed_fq_genome.fa $D/seed_fq_reads_0.fq \
    | grep -v '^@' | sort | test_md5 3d94fd98183a54b0a67644f33270c386

nucmer --sam-long /dev/stdout -l 10 <(echo -e ">101\nggtttatgcgctgttatgtctatggacaaaaaggctacgagaaactgtagccccgttcgctcggacccgcgtcattcgtcggcccagctctacccg") <(echo -e ">21\nggtttatgcgctgttttgtctatggaaaaaaggctacgagaaactgtagccccgttcgctcggtacccgcgtcattcgtcggcccatctctacccg") \
    | grep -v '^@' | test_md5 f656b26b59de04e7c94c7c0c0f7e3a0c
