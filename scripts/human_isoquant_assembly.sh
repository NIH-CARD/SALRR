#!/usr/bin/env bash

usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

            # Your options go here.
            --infile) assert_argument "$1" "$opt"; INFILE="$1"; shift;;
            --outdir) assert_argument "$1" "$opt"; OUTDIR="$1"; shift;;
            --genome) assert_argument "$1" "$opt"; GENOME="$1"; shift;;
            --genedb) assert_argument "$1" "$opt"; GENEDB="$1"; shift;;
            --prefix) assert_argument "$1" "$opt"; PREFIX="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
            --mode) assert_argument "$1" "$opt"; MODE="$1"; shift;;
      
            # Arguments processing. You may remove any unneeded line after the 1st.
            -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
            -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
            *)          usage_error "this should NEVER happen ($opt)";;               # sanity test for previous patterns
    
        esac
    done
    shift  # $EOL
fi



# Rest of code

# making the output directory if it does not exist already

mkdir -p "${OUTDIR}"

# IsoQuant quantification or discovery mode optionality
ISOQUANT_OPTS="--check_canonical --sqanti_output"
if [ "$MODE" = "quantification" ]; then
    ISOQUANT_OPTS="--no_model_construction"  # Quantification only, no novel transcripts
fi

# Run Isoquant 
isoquant.py \
    -t "${THREADS}" \
    --reference ${GENOME} \
    --genedb ${GENEDB} \
    --complete_genedb \
    --bam ${INFILE} \
    --data_type nanopore \
    --prefix ${PREFIX} \
    --count_exons \
    ${ISOQUANT_OPTS} \
    --output ${OUTDIR}
