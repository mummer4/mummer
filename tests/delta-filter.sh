time nucmer --maxmatch --delta ori.delta $D/small_reads_1.fa $D/small_reads_0.fa

# No filter -> pass through
cmp ori.delta <(delta-filter ori.delta)

# Some filtering: check coords
delta-filter -i 95 ori.delta > filter_i95.delta
cmp <(show-coords -H -T -I 95 ori.delta) <(show-coords -H -T filter_i95.delta)
delta-filter -l 95 ori.delta > filter_l95.delta
cmp <(show-coords -H -T -L 95 ori.delta) <(show-coords -H -T filter_l95.delta)
