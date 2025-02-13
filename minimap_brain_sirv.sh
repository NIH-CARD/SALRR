#!/usr/bin/env bash

#SBATCH --cpus-per-task=120
#SBATCH --mem=400g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:100

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify path to sample sheet
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/data/CARD_AUX/LRS_temp/NABEC_RNA/SEQ_REPORTS/basecalling_sample_sheet_dec_24.txt'
# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

basecalling_model="cdna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory paths

BASE_DIR="/data/CARD_AUX/LRS_temp/NABEC_RNA"
MAPPED_DIR="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
PYCHOPPER_DIR="${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}"
ONT_UBAM_DIR="${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}"
SIRV_REF_DIR="/data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV"
REF_DIR="/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
# make output directories and parent if necessary

mkdir -p ${BASE_DIR}/MAPPED/${SAMPLE_ID}/

# load modules

module load samtools/1.21
module load minimap2/2.28

# Output Bam files
SIRV_UNFILTERED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam"
SIRV_FILTERED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam"
UNMAPPED_BAM="${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam"
UNMAPPED_FASTQ="${PYCHOPPER_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.fastq"
BRAIN_MAPPED_BAM="${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam"


# Mapping
# debugging output with output path for unmapped BAM
echo "${SIRV_UNFILTERED_BAM}"
echo "${SIRV_FILTERED_BAM}"
echo "${UNMAPPED_BAM}"
echo "${UNMAPPED_FASTQ}"
echo "${BRAIN_MAPPED_BAM}"
minimap2 \
-t $SLURM_CPUS_PER_TASK \
-ax splice \
--splice-flank=no \
"${SIRV_REF_DIR}/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta" \
${PYCHOPPER_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.trimmed.fastq - \
| samtools view -b - \
| samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> "${SIRV_UNFILTERED_BAM}"
samtools index "${SIRV_UNFILTERED_BAM}" && \
samtools view \
-q 40 \
-F 2304 \
-b \
"${SIRV_UNFILTERED_BAM}" | samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> "${SIRV_FILTERED_BAM}"
samtools index "${SIRV_FILTERED_BAM}" && \
samtools view \
-f 4 \
-b \
"${SIRV_UNFILTERED_BAM}" | samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> "${UNMAPPED_BAM}"
samtools index "${UNMAPPED_BAM}" && \
samtools fastq \
-T* \
-@ $SLURM_CPUS_PER_TASK \
-n \
"${UNMAPPED_BAM}" \
> "${UNMAPPED_FASTQ}" && \
minimap2 \
-t $SLURM_CPUS_PER_TASK \
-ax splice \
-k 14 \
-uf \
{$REF_DIR} \
"${UNMAPPED_FASTQ}" - \
| samtools view \
-q 40 \
-F 2304 \
-b - \
| samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> "${BRAIN_MAPPED_BAM}"
samtools index \
"${BRAIN_MAPPED_BAM}"
