#!/usr/bin/env bash

#SBATCH --cpus-per-task=60
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=04:00:00
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
STR_BRAIN="/data/CARDPB/data/LRS_RNA/data/NGD"
STR_BRAIN="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
PYCHOPPER_DIR="${BASE_DIR}/PYCHOPPER/${SAMPLE_ID}"
BRAIN_UBAM_DIR="${BASE_DIR}/BRAIN_UBAM/${SAMPLE_ID}" # NGD directory for unmapped brain reads
ONT_UBAM_DIR="${BASE_DIR}/ONT_UBAM/${FLOWCELL}" # NGD directory for unmapped reads
REF_GTF="/data/CARDPB/resources/hg38/gencode.v43.annotation.gtf"
REF_FASTA="/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# make output directories and parent if necessary

mkdir -p ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/ #stringtie output directory

mkdir -p ${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/ #isoquant output directory

# load modules
module load stringtie/2.2.3
module load isoquant/3.6.2

# Assembly
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf"
echo "${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_ISO_asm.gtf"
stringtie \
${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
-L \
-p $SLURM_CPUS_PER_TASK \
-G /data/CARDPB/resources/hg38/gencode.v43.annotation.gtf \
-o ${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf \
&& \
isoquant.py \
-t 60 \
--reference ${REF_FASTA} \
--genedb /data/CARDPB/resources/hg38/gencode.v43.annotation.gtf \
--complete_genedb \
--bam ${BASE_DIR}/MAPPED/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
--data_type nanopore \
--check_canonical \
--sqanti_output \
--prefix ${SAMPLE_ID} \
-o ${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}/
# How to run the script
# s b a t c h --array=1-10 script_name.sh (this is an example if you have 10 samples)

