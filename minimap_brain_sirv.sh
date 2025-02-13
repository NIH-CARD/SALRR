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
#
module load samtools/1.21
module load minimap2/2.28

# Mapping
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam"
echo "${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam"
echo "${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam"
echo "${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.fastq"
echo "${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam"
minimap2 \
-t $SLURM_CPUS_PER_TASK \
-ax splice \
--splice-flank=no /data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV/SIRV_ERCC_longSIRV_multi-fasta
_20210507.fasta ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.t
rimmed.fastq - \
| samtools view -b - \
| samtools sort -@ $SLURM_CPUS_PER_TASK - \
> ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam
samtools index ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam && \
samtools view \
-q 40 \
-F 2304 \
-b \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam | samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam
samtools index ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam && \
samtools view \
-f 4 \
-b \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRVome_mapped_unfiltered.sorted.bam | samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam
samtools index ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam && \
samtools fastq \
-T* \
-@ $SLURM_CPUS_PER_TASK \
-n \
${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.sorted.bam \
> ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_unmapped.fastq && \
minimap2 \
-t $SLURM_CPUS_PER_TASK \
-ax splice \
-k 14 \
-uf \
/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/${S
AMPLE_ID}_${FLOWCELL}_brain_unmapped.fastq - \
| samtools view \
-q 40 \
-F 2304 \
-b - \
| samtools sort \
-@ $SLURM_CPUS_PER_TASK - \
> ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam
samtools index \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam
