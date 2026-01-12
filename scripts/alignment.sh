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
            --genome) assert_argument "$1" "$opt"; GENOME="$1"; shift;;
            --sirvome) assert_argument "$1" "$opt"; SIRV_REF="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
            --skip-sirv) SKIP_SIRV=true; shift;;
      
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

# Default: SIRV analysis enabled unless --skip-sirv flag is set
SKIP_SIRV=${SKIP_SIRV:-false}

# making the output directory if it does not exist

mkdir -p "${OUTDIR}"

if [ "$SKIP_SIRV" = "false" ]; then
    echo "Running with SIRV spike-in analysis..."
    
    # map reads to SIRVome
    minimap2 \
        -t "${THREADS}" \
        -ax splice \
        --splice-flank=no \
        "${SIRV_REF}" \
        "${INFILE}" - \
        | samtools view -b - \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam"
    
    # index bam, exit if it fails
    samtools index "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" || exit 1

    # filter the SIRV mapped reads
    samtools view \
        -q 40 \
        -F 2304 \
        -b \
        "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam"

    # index the mapped filtered sirv bam reads, exit if it fails
    samtools index "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam" || exit 1

    # convert unmapped human reads to fastq, exit if it fails
    samtools view \
        -f 4 \
        -b \
        "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%_human_mapped.sorted.bam}_human_unmapped.sorted.bam" || exit 1

    # index the bam, exit if it fails
    samtools index "${OUTFILE%_human_mapped.sorted.bam}_human_unmapped.sorted.bam" || exit 1

    # convert unmapped human reads to fastq, exit if it fails
    samtools fastq \
        -T* \
        -@ "${THREADS}" \
        -n \
        "${OUTFILE%_human_mapped.sorted.bam}_human_unmapped.sorted.bam" \
        > "${INFILE%.trimmed.fastq}_human_unmapped.fastq" || exit 1

    # map the human reads to the genome
    minimap2 \
        -t "${THREADS}" \
        -ax splice \
        "${GENOME}" \
        "${INFILE%.trimmed.fastq}_human_unmapped.fastq" - \
        | samtools view \
        -q 40 \
        -F 2304 \
        -b - \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE}" || exit 1

    # index the human mapped bam, exit if it fails
    samtools index "${OUTFILE}" || exit 1
    
    # # Generate mapping statistics with SIRV
    # STATS_FILE="${OUTDIR}/$(basename ${OUTFILE%.bam})_mapping_stats.txt"
    # TOTAL_READS=$(echo $(cat "${INFILE}" | wc -l) / 4 | bc)
    # SIRV_MAPPED=$(samtools view -c -F 4 "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam")
    # HUMAN_MAPPED=$(samtools view -c "${OUTFILE}")
    # SIRV_UNMAPPED=$(samtools view -c "${OUTFILE%_human_mapped.sorted.bam}_human_unmapped.sorted.bam")
    # UNMAPPED=$((TOTAL_READS - SIRV_MAPPED - HUMAN_MAPPED))
    
    # SIRV_PCT=$(echo "scale=2; $SIRV_MAPPED * 100 / $TOTAL_READS" | bc)
    # HUMAN_PCT=$(echo "scale=2; $HUMAN_MAPPED * 100 / $TOTAL_READS" | bc)
    # UNMAPPED_PCT=$(echo "scale=2; $UNMAPPED * 100 / $TOTAL_READS" | bc)
    # HUMAN_OF_NONSIRVS=$(echo "scale=2; $HUMAN_MAPPED * 100 / $SIRV_UNMAPPED" | bc)
    
    # echo -e "Total Reads:\t$TOTAL_READS" > "$STATS_FILE"
    # echo -e "SIRV Mapped:\t$SIRV_MAPPED\t($SIRV_PCT%)" >> "$STATS_FILE"
    # echo -e "Human Mapped:\t$HUMAN_MAPPED\t($HUMAN_PCT%)" >> "$STATS_FILE"
    # echo -e "Unmapped:\t$UNMAPPED\t($UNMAPPED_PCT%)" >> "$STATS_FILE"
    # echo -e "Human % of non-SIRV reads:\t$HUMAN_OF_NONSIRVS%" >> "$STATS_FILE"
    
else
    echo "Skipping SIRV analysis, mapping directly to human genome..."
    
    # Map directly to human genome
    minimap2 \
        -t "${THREADS}" \
        -ax splice \
        "${GENOME}" \
        "${INFILE}" - \
        | samtools view \
        -q 40 \
        -F 2304 \
        -b - \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE}" || exit 1
    
    # index the human mapped bam, exit if it fails
    samtools index "${OUTFILE}" || exit 1
    
    # # Generate mapping statistics without SIRV
    # STATS_FILE="${OUTDIR}/$(basename ${OUTFILE%.bam})_mapping_stats.txt"
    # TOTAL_READS=$(echo $(cat "${INFILE}" | wc -l) / 4 | bc)
    # HUMAN_MAPPED=$(samtools view -c "${OUTFILE}")
    # UNMAPPED=$((TOTAL_READS - HUMAN_MAPPED))
    
    # HUMAN_PCT=$(echo "scale=2; $HUMAN_MAPPED * 100 / $TOTAL_READS" | bc)
    # UNMAPPED_PCT=$(echo "scale=2; $UNMAPPED * 100 / $TOTAL_READS" | bc)
    
    # echo -e "Total Reads:\t$TOTAL_READS" > "$STATS_FILE"
    # echo -e "Human Mapped:\t$HUMAN_MAPPED\t($HUMAN_PCT%)" >> "$STATS_FILE"
    # echo -e "Unmapped:\t$UNMAPPED\t($UNMAPPED_PCT%)" >> "$STATS_FILE"
    # echo -e "SIRV analysis: DISABLED" >> "$STATS_FILE"
fi