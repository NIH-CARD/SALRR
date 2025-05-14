#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=10g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=04:00:00
#SBATCH --gres=lscratch:50
#SBATCH --partition=quick

module load python

# Define base directory
BASE_DIR="/data/CARDPB/data/LRS_RNA/data/NGD"

# make output directories and parent if necessary
mkdir -p ${BASE_DIR}/SIRV_analysis/

# Sample sheet path
SAMPLE_SHEET="${BASE_DIR}/SAMPLE_SHEETS/NDG_sample_sheet.tsv"

# Get sample ID and flowcell from sample sheet
N=${SLURM_ARRAY_TASK_ID}
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# Construct GTF file path
GTF_FILE="${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf"

# Debugging input gtf file
echo "Processing GTF file at: $GTF_FILE"

# Run the Python script
python sirvsuite_gtf_parsing_main.py "$GTF_FILE" "$SAMPLE_ID" "$FLOWCELL"