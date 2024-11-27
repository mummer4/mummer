#! /usr/bin/env bash

set -e

libexec="@LIBEXEC_DIR@"

usage() {
    echo "Usage: $0 [-h,--help] command [options] [args....]"
}

help() {
    cat <<EOF
Run a tool from the MUMmer toolbox.

Subcommands:
mummer       : find MUMs and MEMs between sequences
nucmer       : nucleotide aligner
promer       : protein aligner
mummerplot   : plot 2D alignments
show-coords  : show the coordinates of alignments
show-aligns  : show detailed alignment information
show-diff    : show diff
show-snps    : show snps
show-tiling  : show tiling
delta-filter : filter a delta file
EOF
}

# If ask for global help
if [ "$1" == "-h" -o "$1" == "--help" ]; then
    usage
    echo
    help
    exit 0
fi

command=$(basename "$0")

case "$command" in
    (*mumtool*) # Called as 'mumtool subcommand ....'
        if [ "$#@" == "1" ]; then
            { echo "Missing subcommand"; usage; } >&2
            exit 1
        fi
        command=$1
        shift 1
        ;;

    (*) # Called as a symlink with the name of the subcommand
        ;;
esac

program="${libexec}/${command}"
if ! [ -x "$program" ]; then
    { echo "Invalid subcommand '$command'" >&2; usage; } >&2
    exit 1
fi

exec "$program" "$@"
