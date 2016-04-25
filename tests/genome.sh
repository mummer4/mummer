time  sed -e 's/^>\([^[:space:]]\+\).*/>\1/' $D/seed_reads_2.fa | tee genome | nucmer -G --delta /dev/stdout $D/seed_genome.fa /dev/stdin | \
    tee genome.delta | tail -n +3 | test_md5 8328b1577d8656eaa53aa61a113d89b0
