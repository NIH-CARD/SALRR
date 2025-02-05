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

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory paths
#BASE_DIR=/data/CARD_AUX/LRS_temp/December_24/ # for example
BASE_DIR="/data/CARD_AUX/LRS_temp/NABEC_RNA"

# make output directories and parent if necessary
mkdir -p ${BASE_DIR}/U_FASTQ/${SAMPLE_ID}
mkdir -p ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID} #for pychopper output
mkdir -p ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/stats #for pychopper summary stats

# load modules
#
module load samtools/1.21
module load pychopper/2.7.10

# Trimming
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/U_FASTQ/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.trimmed.fastq"
samtools fastq -T* -@ $SLURM_CPUS_PER_TASK -n ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.bam > ${BASE_DIR}/U_FASTQ/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.fastq && \
pychopper -t $SLURM_CPUS_PER_TASK -m phmm -k PCS114 -r ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/stats/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.pdf -S ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/stats/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.tsv ${BASE_DIR}/U_FASTQ/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.fastq ${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_cdna_r10.4.1_e8.2_400bps_sup@v5.0.0.trimmed.fastq 


#How to run the script
#s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)
