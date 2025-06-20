#!/usr/bin/env bash

#SBATCH --partition=gpu
#SBATCH --cpus-per-task=30
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:50,gpu:a100:2

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

#flow cell specific basecalling model
basecalling_model="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# base directory paths

BASE_DIR=/data/CARD_AUX/LRS_temp/NABEC_RNA/ # for example
ONT_UBAM_DIR="${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}"

# make output directory and parent if necessary
mkdir -p "${ONT_UBAM_DIR}"

# load modules
# change 0.8.1 to 0.9.0 module (not yet default)
module load dorado/0.9.0
module load pod5/0.3.15


# debugging output with output path for unmapped BAM
echo "${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_c${basecalling_model}.bam"

# Basecalling with dorado
dorado basecaller \
    --no-trim \
    --estimate-poly-a \
    -x cuda:all \
    ${DORADO_MODELS}/${basecalling_model} \
    ${BASE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}/${FLOWCELL}/pod5 \
    --skip-model-compatibility-check \
    > ${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_c${basecalling_model}.bam

# debugging output with output path for sequencing summary QC text file
echo "${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_sequencing_summary_v5.0.0.txt"

# Generating sequencing summary reports
dorado summary \
    ${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_c${basecalling_model}.bam \
    > ${ONT_UBAM_DIR}/${SAMPLE_ID}_${FLOWCELL}_sequencing_summary_v5.0.0.txt

# How to run the script
# sbatch --array=1-10 script_name.sh (this is an example if you have 10 samples)
