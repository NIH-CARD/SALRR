#!/usr/bin/env bash
set -euo pipefail

# Re-quantify ONE sample's BAM against the unified cohort GTF with IsoQuant.
# Output nests as <quant_dir>/<sample_id>/<sample_id>.<type>.tsv

usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

              # Command line options below
              # Single values
              --bam)          assert_argument "$1" "$opt"; BAM="$1"; shift ;;
              --sample_id)    assert_argument "$1" "$opt"; SAMPLE_ID="$1"; shift ;;
              --cohort_gtf)   assert_argument "$1" "$opt"; COHORT_GTF="$1"; shift ;;
              --genome)       assert_argument "$1" "$opt"; GENOME="$1"; shift ;;
              --quant_dir)    assert_argument "$1" "$opt"; QUANT_DIR="$1"; shift ;;
              --threads)      assert_argument "$1" "$opt"; THREADS="$1"; shift ;;

              # Arguments processing. You may remove any unneeded line after the 1st.
              -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
              --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
              -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
              --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
              -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
              *)          usage_error "this should NEVER happen ($opt)";;

        esac
    done
    shift   # $EOL
fi

# Rest of code

mkdir -p "${QUANT_DIR}"

[ -s "${BAM}" ] || { echo "ERROR: missing BAM for ${SAMPLE_ID}: ${BAM} (needed for re-quant)" >&2; exit 1; }

###############################################################################
# IsoQuant re-quantification against the cohort GTF
#   Count this sample's reads against the single unified cohort model so every
#   sample is quantified on the same reference (comparable across the cohort).
###############################################################################

isoquant \
    -t "${THREADS}" \
    --reference "${GENOME}" \
    --genedb "${COHORT_GTF}" \
    --transcript_quantification unique_only \
    --gene_quantification unique_splicing_consistent \
    --no_model_construction \
    --data_type nanopore \
    --count_exons \
    --bam "${BAM}" \
    --prefix "${SAMPLE_ID}" \
    -o "${QUANT_DIR}"