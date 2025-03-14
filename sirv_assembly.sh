#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=30g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=04:00:00
#SBATCH --gres=lscratch:50
#SBATCH --partition=quick

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
STR_ASM="${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV"
ISO_ASM="${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/SIRV"
MAPPED_DIR="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
ref_base_dir="/data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV"
REF_GTF="${ref_base_dir}/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf"
CORRECTED_GTF="${ref_base_dir}/SIRV_ERCC_longSIRV_multi-fasta_20210507.corrected.gtf"
REF_FASTA="${ref_base_dir}/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta"


# make output directories and parent if necessary

mkdir -p "${STR_ASM}" "${ISO_ASM}" # stringtie isoquant output directory

# load modules
module load stringtie/2.2.3
module load isoquant/3.6.2

# Assembly
# debugging output with output path for unmapped BAM
echo "${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf"
echo "${ISO_ASM}"
# Run SIRV Stringtie Assembly
stringtie \
    ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
    -L \
    -p $SLURM_CPUS_PER_TASK \
    -G ${REF_GTF} \
    -o ${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf || exit 1

# Run SIRV Isoquant Assembly

isoquant.py \
    -t 60 \
    --reference ${REF_FASTA} \
    --genedb ${CORRECTED_GTF} \
    --bam ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
    --data_type nanopore \
    --check_canonical \
    --sqanti_output \
    --prefix ${SAMPLE_ID} \
    -o ${ISO_ASM} || exit 1

# Run SIRV Stringtie quantification

stringtie \
    ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
    -e \
    -L \
    -p $SLURM_CPUS_PER_TASK \
    -G ${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf \
    -o ${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_sirv_quantified_STR_asm.gtf


# How to run the script
# s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)

