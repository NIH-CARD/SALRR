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
            --qc-dir) assert_argument "$1" "$opt"; QC_DIR="$1"; shift;;
            --qc-file) assert_argument "$1" "$opt"; QC_FILE="$1"; shift;;
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

# making a directory for output QC if it does not exist

mkdir -p "${QC_DIR}"

if [ "$SKIP_SIRV" = "false" ]; then
    echo "Running with SIRV spike-in mapping..."
    
 ############# MAPPING WITH SIRV SPIKE-IN #############


    # map reads to SIRVome (unfiltered)
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

    #Generating QC stats of SIRV mapping
    samtools stats \
        -@ "${THREADS}" \
        "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" \
        > "${QC_FILE%_stats.txt}_sirv_stats.txt" || exit 1

    # filter the SIRV mapped reads for downstream SIRV analysis
    samtools view \
        -q 40 \
        -F 2304 \
        -b \
        "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam"

    # index the mapped filtered sirv bam reads, exit if it fails
    samtools index "${OUTFILE%_human_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam" || exit 1

    # Extract unmapped reads from SIRV (these should be human or sample only reads)
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


    # Map human reads (SIRV-unmapped) to human genome (unfiltered for accurate QC stats)
    minimap2 \
        -t "${THREADS}" \
        -ax splice \
        "${GENOME}" \
        "${INFILE%.trimmed.fastq}_human_unmapped.fastq" - \
        | samtools view \
        -b - \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam" || exit 1

    # Index the human unfiltered mapped bam, exit if it fails
    samtools index "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam" || exit 1

    # Human QC stats (on unfiltered SIRV-depleted - shows % of SIRV-depleted reads that map to human)
    samtools stats \
        -@ "${THREADS}" \
        "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam" \
        > "${QC_FILE%_stats.txt}_human_stats.txt" || exit 1

     # Filter human mapped reads (final BAM for assembly/quantification)
    samtools view \
        -q 40 \
        -F 2304 \
        -b \
        "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam" \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE}" || exit 1

    # Index the final human mapped filtered bam, exit if it fails
    samtools index "${OUTFILE}" || exit 1


    # Final BAM stats (shows impact of quality filtering)
    samtools stats \
        -@ "${THREADS}" \
        "${OUTFILE}" \
        > "${QC_FILE%_stats.txt}_human_filtered_stats.txt" || exit 1

    # Cleaning up the unfiltered BAM to clear disk space (these files are too large to keep)
    echo "Removing intermediate files..."
    rm "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam"
    rm "${OUTFILE%_human_mapped.sorted.bam}_human_unfiltered.sorted.bam.bai"
    
    echo "Alignment with SIRV spike-in complete!"



 ############# MAPPING WITH NO SIRV ##################

else
    echo "Skipping SIRV alignment, mapping directly to human genome..."
    

    # Map directly to human genome
    minimap2 \
        -t "${THREADS}" \
        -ax splice \
        "${GENOME}" \
        "${INFILE}" - \
        | samtools view \
        -b - \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam" || exit 1
    
    # index the human mapped bam, exit if it fails
    samtools index "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam" || exit 1

    # Human QC stats (on unfiltered)
    samtools stats \
        -@ "${THREADS}" \
        "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam" \
        > "${QC_FILE%_stats.txt}_human_stats.txt" || exit 1
    
    # Filter human mapped reads (final BAM)
    samtools view \
        -q 40 \
        -F 2304 \
        -b \
        "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam" \
        | samtools sort \
        -@ "${THREADS}" - \
        > "${OUTFILE}" || exit 1
    
    # Index the final human mapped filtered bam, exit if it fails
    samtools index "${OUTFILE}" || exit 1

    # Final BAM stats
    samtools stats \
        -@ "${THREADS}" \
        "${OUTFILE}" \
        > "${QC_FILE%_stats.txt}_human_filtered_stats.txt" || exit 1


    #  Cleaning up the unfiltered BAM to clear disk space (these files are too large to keep)
    echo "Removing intermediate files..."
    rm "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam"
    rm "${OUTFILE%.sorted.bam}_unfiltered.sorted.bam.bai"
    
    echo "Human mapping complete!"



fi