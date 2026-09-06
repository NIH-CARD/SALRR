#!/usr/bin/env bash
set -euo pipefail

# Build one unified cohort GTF from per-sample annotated GTF files with IsoMatch.

usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

              # Command line options below
              # matrix table 3-column tsv sample id | annotated gtf path | bam path all tab separated
              --manifest)       assert_argument "$1" "$opt"; MANIFEST="$1"; shift ;; 
              # Single values
              --cohort_dir)     assert_argument "$1" "$opt"; COHORT_DIR="$1"; shift ;;
              --cohort_name)    assert_argument "$1" "$opt"; COHORT_NAME="$1"; shift ;;
              --genome)         assert_argument "$1" "$opt"; GENOME="$1"; shift ;;
              --ref_gtf)        assert_argument "$1" "$opt"; REF_GTF="$1"; shift ;;

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

# Derived layout
COHORT_GTF_DIR="${COHORT_DIR}/gtf"
COHORT_WORK_DIR="${COHORT_DIR}/work"
COHORT_GTF="${COHORT_GTF_DIR}/${COHORT_NAME}.annotated.gtf"

mkdir -p "${COHORT_GTF_DIR}" "${COHORT_WORK_DIR}"

###############################################################################
# isomatch cross-sample merge
#   Sort + index each per-sample annotated GTF, merge them all into one,
#   then classify the merged model against the reference annotation.
###############################################################################

# Parse the 3-column tab-separated table line by line to get sample ID, annotated GTF path, and BAM path
SORTED_GTFS=()
while IFS=$'\t' read -r sample_id annotated_gtf human_mapped_bam; do
    [ -z "${sample_id}" ] && continue          # skip blank lines
    sorted_gtf="${COHORT_WORK_DIR}/${sample_id}.sorted.gtf"

    bedtools sort -i "${annotated_gtf}" > "${sorted_gtf}"
    isomatch index --ref-fa "${GENOME}" "${sorted_gtf}"
    SORTED_GTFS+=("${sorted_gtf}")
done < "${MANIFEST}"

# Merge all per-sample GTFs into a single cohort model
COHORT_MERGE_PREFIX="${COHORT_GTF_DIR}/${COHORT_NAME}_merged"
isomatch merge \
    --ref-fa "${GENOME}" \
    -o "${COHORT_MERGE_PREFIX}" \
    "${SORTED_GTFS[@]}"

# Classify (annotate) the merged model against the reference
COHORT_CLASSIFY_PREFIX="${COHORT_GTF_DIR}/${COHORT_NAME}"
isomatch index --ref-fa "${GENOME}" "${COHORT_MERGE_PREFIX}.merged.gtf.gz"
isomatch classify \
    --ref-gtf "${REF_GTF}" \
    --ref-fa "${GENOME}" \
    --out "${COHORT_CLASSIFY_PREFIX}" \
    "${COHORT_MERGE_PREFIX}.merged.gtf.gz"

# Decompress the annotated cohort GTF to the final plain-text deliverable.
zcat "${COHORT_CLASSIFY_PREFIX}.annotated.gtf.gz" > "${COHORT_GTF}"
