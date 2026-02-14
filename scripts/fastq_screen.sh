#!/usr/bin/env bash

#Argument parsing
usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }

# # Default: screen a subset of 100000 reads (set to 0 for all reads)
# SUBSET=100000

if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in
            --infile)  assert_argument "$1" "$opt"; INFILE="$1"; shift;;
            --outdir)  assert_argument "$1" "$opt"; OUTDIR="$1"; shift;;
            --conf)    assert_argument "$1" "$opt"; CONF="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
            --subset)  assert_argument "$1" "$opt"; SUBSET="$1"; shift;;

            -|''|[!-]*) set -- "$@" "$opt";;
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;
            -*)         usage_error "unknown option: '$opt'";;
            *)          usage_error "this should NEVER happen ($opt)";;
        esac
    done
    shift  # $EOL
fi

# Validate required arguments
if [[ -z "${INFILE}" || -z "${OUTDIR}" || -z "${CONF}" ]]; then
    echo "ERROR: --infile, --outdir, and --conf are required"
    exit 1
fi

mkdir -p "${OUTDIR}"

echo "=== FastQ Screen ==="
echo "Input:   ${INFILE}"
echo "Output:  ${OUTDIR}"
echo "Config:  ${CONF}"
echo "Threads: ${THREADS}"
echo "Subset:  ${SUBSET} (0 = all reads)"

# Run fastq_screen with minimap2 aligner
fastq_screen \
    --aligner minimap2 \
    --conf "${CONF}" \
    --threads "${THREADS}" \
    --subset "${SUBSET}" \
    --outdir "${OUTDIR}" \
    "${INFILE}"

echo "=== FastQ Screen complete ==="
