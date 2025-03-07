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
#BASE_DIR=/data/CARD_AUX/LRS_temp/December_24/ # for example
BASE_DIR="/data/CARDPB/data/LRS_RNA/data/NGD"

# make output directories and parent if necessary

mkdir -p ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/ #stringtie output directory

mkdir -p ${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/SIRV/ #isoquant output directory

# load modules
module load stringtie/2.2.3
module load isoquant/3.6.2

# Assembly
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf"
echo "${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/SIRV/"
stringtie \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
-L \
-p $SLURM_CPUS_PER_TASK \
-G /data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf \
-o ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf \
&& \
isoquant.py \
-t 60 \
--reference /data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta \
--genedb /data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/REFERENCES/SIRV/SIRV_ERCC_longSIRV_multi-fasta_20210507.corrected.gtf \
--bam ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
--data_type nanopore \
--check_canonical \
--sqanti_output \
--prefix ${SAMPLE_ID} \
-o ${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/SIRV/ \
&& \
stringtie \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_SIRV_mapped_filtered.sorted.bam \
-e \
-L \
-p $SLURM_CPUS_PER_TASK \
-G ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/${SAMPLE_ID}_${FLOWCELL}_sirv_STR_asm.gtf \
-o ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/SIRV/${SAMPLE_ID}_${FLOWCELL}_sirv_quantified_STR_asm.gtf


# How to run the script
# s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)

