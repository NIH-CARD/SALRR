#!/usr/bin/env bash

#SBATCH --cpus-per-task=120
#SBATCH --mem=200g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:100
#SBATCH --job-name=minimap_brain_sirv
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
# flow cell specific basecalling model
basecalling_model="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory paths

BASE_DIR="/data/CARDPB/data/LRS_RNA/data/NGD"
MAPPED_DIR="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
PYCHOPPER_DIR="${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}"
BRAIN_UBAM_DIR="${BASE_DIR}/BRAIN_UBAM/${SAMPLE_ID}" # NGD directory for unmapped brain reads
ONT_UBAM_DIR="${BASE_DIR}/ONT_UBAM/${FLOWCELL}" # NGD directory for unmapped reads
SIRV_REF_DIR="/data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV"
REF_DIR="/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
# make output directories and parent if necessary

mkdir -p "${MAPPED_DIR}" "${BRAIN_UBAM_DIR}"

# load modules

module load samtools/1.21
module load minimap2/2.28

# Output Bam files
SIRV_UNFILTERED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam"
SIRV_FILTERED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam"
UNMAPPED_BAM="${BRAIN_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam"
UNMAPPED_FASTQ="${PYCHOPPER_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.fastq"
BRAIN_MAPPED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam"

# debugging output with output path for unmapped BAM
echo "${SIRV_UNFILTERED_BAM}"
echo "${SIRV_FILTERED_BAM}"
echo "${UNMAPPED_BAM}"
echo "${UNMAPPED_FASTQ}"
echo "${BRAIN_MAPPED_BAM}"

# Mapping
# 1. Map reads to SIRVome
# 2. Filter SIRV mapped reads
# 3. Extract brain unmapped reads
# 4. Map brain reads to genome

# map reads to SIRVome
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -ax splice \
    --splice-flank=no \
    "${SIRV_REF_DIR}/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta" \
    ${PYCHOPPER_DIR}/${SAMPLE_ID}_${FLOWCELL}_c${basecalling_model}.trimmed.fastq - \
    | samtools view -b - \
    | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${SIRV_UNFILTERED_BAM}"

# index the bam, exit if it fails
samtools index "${SIRV_UNFILTERED_BAM}" || exit 1

# filter the SIRV mapped reads
samtools view \
    -q 40 \
    -F 2304 \
    -b \
    "${SIRV_UNFILTERED_BAM}" | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${SIRV_FILTERED_BAM}"

# index the bam, exit if it fail   
samtools index "${SIRV_FILTERED_BAM}" || exit 1

# extract unmapped brain reads, exit if it fails
samtools view \
    -f 4 \
    -b \
    "${SIRV_UNFILTERED_BAM}" | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${UNMAPPED_BAM}" || exit 1

# index the bam, exit if it fails
samtools index "${UNMAPPED_BAM}" || exit 1

# convert unmapped brain reads to fastq, exit if it fails
samtools fastq \
    -T* \
    -@ $SLURM_CPUS_PER_TASK \
    -n \
    "${UNMAPPED_BAM}" \
    > "${UNMAPPED_FASTQ}" || exit 1

# map the brain reads to the genome
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -ax splice \
    -k 14 \
    -uf \
    "${REF_DIR}" \
    "${UNMAPPED_FASTQ}" - \
    | samtools view \
    -q 40 \
    -F 2304 \
    -b - \
    | samtools sort \
    -@ $SLURM_CPUS_PER_TASK - \
    > "${BRAIN_MAPPED_BAM}" || exit 1

# index the brain mapped bam, exit if it fails
samtools index \
    "${BRAIN_MAPPED_BAM}" || exit 1