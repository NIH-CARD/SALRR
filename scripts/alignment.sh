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

# loading modules
module load minimap2/2.29
module load samtools/1.21

# making the output directory if it does not exist

mkdir -p "${OUTDIR}"

# map reads to SIRVome
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -ax splice \
    --splice-flank=no \
    "${SIRV_REF}" \
    "${INFILE}" - \
    | samtools view -b - \
    | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam"

# index the bam, exit if it fails
samtools index "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" || exit 1

# filter the SIRV mapped reads
samtools view \
    -q 40 \
    -F 2304 \
    -b \
    "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam"

# index the mapped filtered sirv bams, exit if it fails
samtools index "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_filtered.sorted.bam" || exit 1

# extract unmapped brain reads, exit if it fails
samtools view \
    -f 4 \
    -b \
    "${OUTFILE%_brain_mapped.sorted.bam}_SIRVome_mapped_unfiltered.sorted.bam" | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${OUTFILE%_brain_mapped.sorted.bam}_brain_unmapped.sorted.bam" || exit 1

# index the bam, exit if it fails
samtools index "${OUTFILE%_brain_mapped.sorted.bam}_brain_unmapped.sorted.bam" || exit 1

# convert unmapped brain reads to fastq, exit if it fails
samtools fastq \
    -T* \
    -@ $SLURM_CPUS_PER_TASK \
    -n \
    "${OUTFILE%_brain_mapped.sorted.bam}_brain_unmapped.sorted.bam" \
    > "${INFILE%.trimmed.fastq}_brain_unmapped.fastq" || exit 1

# map the brain reads to the genome
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -ax splice \
    -k 14 \
    -uf \
    "${GENOME}" \
    "${INFILE%.trimmed.fastq}_brain_unmapped.fastq" - \
    | samtools view \
    -q 40 \
    -F 2304 \
    -b - \
    | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${OUTFILE}" || exit 1

# index the brain mapped bam, exit if it fails
samtools index \
    "${OUTFILE}" || exit 1
