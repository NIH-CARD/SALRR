#!/usr/bin/env bash




# Argument parsing
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
            --outfile) assert_argument "$1" "$opt"; OUTFILE="$1"; shift;;
            --outdir) assert_argument "$1" "$opt"; OUTDIR="$1"; shift;;
            --qc-file) assert_argument "$1" "$opt"; QC_FILE="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
            --model) assert_argument "$1" "$opt"; MODEL="$1"; shift;;
      
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

# Deriving QC directory from QC file path
QC_DIR=$(dirname "${QC_FILE}")

# Create output directory if it does not exist
mkdir -p "${OUTDIR}" "${QC_DIR}"


# Basecalling with dorado
dorado basecaller \
    --no-trim \
    --estimate-poly-a \
    -x cuda:all \
    ${DORADO_MODELS}/${MODEL} \
    ${INFILE} \
    --skip-model-compatibility-check \
    > ${OUTFILE} || exit 1


# Generate sequencing summary reports
dorado summary \
    ${OUTFILE} \
    > ${OUTFILE%.bam}.summary.txt || exit 1


# Generate QC stats after basecalling prior to mapping
cramino \
    ${OUTFILE} \
    --threads ${THREADS} \
    --ubam \
    --checksum \
    --spliced \
    --hist \
    > ${QC_FILE} || exit 1