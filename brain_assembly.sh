#!/usr/bin/env bash

#SBATCH --cpus-per-task=60
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:50
#SBATCH --partition=norm

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify path to sample sheet
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/data/CARDPB/data/LRS_RNA/data/NGD/SAMPLE_SHEETS/NDG_sample_sheet.tsv'
# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory paths
BASE_DIR="/data/CARDPB/data/LRS_RNA/data/NGD"
STR_ASM="${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}"
ISO_ASM="${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}"
MAPPED_DIR="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
REF_GTF="/data/CARDPB/resources/hg38/gencode.v43.annotation.gtf"
REF_FASTA="/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# make output directories and parent if necessary

mkdir -p "${STR_ASM}" "${ISO_ASM}" # Stringtie and IsoQuant output directory

# load modules
module load stringtie/2.2.3
module load isoquant/3.6.2

# Assembly
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf"
echo "${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_ISO_asm.gtf"

# Stringtie Assembly
stringtie \
    ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
    -L \
    -p $SLURM_CPUS_PER_TASK \
    -G ${REF_GTF} \
    -o ${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf || exit 1

# Isoquant Assembly

isoquant.py \
    -t 60 \
    --reference ${REF_FASTA} \
    --genedb ${REF_GTF} \
    --complete_genedb \
    --bam ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
    --data_type nanopore \
    --check_canonical \
    --sqanti_output \
    --prefix ${SAMPLE_ID} \
    --count_exons \
    -o ${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID} || exit 1
# How to run the script
# s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)

