#!/usr/bin/env bash

#SBATCH --cpus-per-task=60
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:50

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
ONT_BAM_DIR="${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}"
FASTQ_DIR="${BASE_DIR}/U_FASTQ/${SAMPLE_ID}"
PYCHOPPER_DIR="${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}"
STATS_DIR="${PYCHOPPER_DIR}/stats"


# make output directories and parent if necessary

mkdir -p "${FASTQ_DIR}" "${PYCHOPPER_DIR}" "${STATS_DIR}"  



# load modules

module load samtools/1.21
module load pychopper/2.7.10

# Defining file paths

INPUT_BAM="${ONT_BAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.bam"
OUTPUT_FASTQ="${FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.fastq"
TRIMMED_FASTQ="${PYCHOPPER_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.trimmed.fastq"
PYCHOPPER_PDF="${STATS_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.pdf"
PYCHOPPER_TSV="${STATS_DIR}/${SAMPLE_ID}_${FLOWCELL}_${basecalling_model}.tsv"

# Debugging output
echo "Processing BAM file: $INPUT_BAM"
echo "Output FASTQ: $OUTPUT_FASTQ"
echo "Trimmed FASTQ: $TRIMMED_FASTQ"

# Converting BAM to FASTQ using samtools

samtools fastq \
-T* \
-@ $SLURM_CPUS_PER_TASK \
-n \
"${INPUT_BAM}" \
> "${OUTPUT_FASTQ}" && \
pychopper \
-t $SLURM_CPUS_PER_TASK \
-m phmm \
-k PCS114 \
-r "${PYCHOPPER_PDF}" \
-S "${PYCHOPPER_TSV}" \
"${OUTPUT_FASTQ}" \
"${TRIMMED_FASTQ}" 


#How to run the script
#s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)
