#!/usr/bin/env bash
set -e




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
            --kit) assert_argument "$1" "$opt"; KIT="$1"; shift;;
            --outdir) assert_argument "$1" "$opt"; OUTDIR="$1"; shift;;
            --qc-dir) assert_argument "$1" "$opt"; QC_DIR="$1"; shift;;
            --qc-file) assert_argument "$1" "$opt"; QC_FILE="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
      
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

#making the output directory and QC directory if they do not exist

mkdir -p "${OUTDIR}" "${QC_DIR}"

# Converting BAM to FASTQ using samtools
samtools fastq \
  -T* \
  -@ "${THREADS}" \
  -n \
  "${INFILE}" \
  > "${OUTFILE%.trimmed.fastq}.fastq" || exit 1

# Running pychopper for trimming

pychopper \
  -t "${THREADS}" \
  -m phmm \
  -k "${KIT}" \
  -r "${OUTFILE%.trimmed.fastq}.pdf" \
  -S "${QC_FILE}" \
  "${OUTFILE%.trimmed.fastq}.fastq" \
  "${OUTFILE}" 


