#!/usr/bin/env bash
#SBATCH --cpus-per-task=60
#SBATCH --mem=300g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=72:00:00
#SBATCH --gres=lscratch:50
#SBATCH --partition=norm
#SBATCH --job-name=tama_merge

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify path to sample sheet
# tab delimited file with list of samples and flow cells in columns 1 and 2
#SAMPLE_SHEET='/data/CARDPB/data/LRS_RNA/data/NGD/SAMPLE_SHEETS/NDG_sample_sheet.tsv'
SAMPLE_SHEET='/data/CARD_AUX/LRS_temp/March_21/SAMPLE_SHEETS/basecalling_sample_sheet_march_21.txt'

# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory paths
BASE_DIR="/data/CARD_AUX/LRS_temp/March_21"
TAMA_BASE_DIR="/data/CARDPB/data/LRS_RNA/projects/cedrics_nabec/tama"
TAMA_CONVERTER="${TAMA_BASE_DIR}/tama_go/format_converter"
STR_ASM="${BASE_DIR}/ASSEMBLY/STRINGTIE/${SAMPLE_ID}"
ISO_ASM="${BASE_DIR}/ASSEMBLY/ISOQUANT/${SAMPLE_ID}"
MERGED_DIR="${BASE_DIR}/ASSEMBLY/MERGED/${SAMPLE_ID}"
MAPPED_DIR="${BASE_DIR}/MAPPED/${SAMPLE_ID}"
REF_GTF="/data/CARDPB/resources/hg38/gencode.v43.annotation.gtf"
REF_FASTA="/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# make output directories and parent if necessary

mkdir -p "${MERGED_DIR}" # Merged assembly output directory

# load modules
module load python/3.10
ml isoquant/3.6.2
#preprocess stringtie gtf to remove the . from strand column
#this is necessary for the merge step to work properly
awk -F '\t' '$7 != "." {print}' ${STR_ASM}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf \
    > ${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf || exit 1

#convert stringtie gtf to bed12
python \
${TAMA_CONVERTER}/tama_format_gtf_to_bed12_stringtie.py \
    ${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.gtf \
    ${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.bed || exit 1

#convert isoquant to bed12
python \
${TAMA_CONVERTER}/tama_format_gtf_to_bed12_ensembl.py \
    ${ISO_ASM}/${SAMPLE_ID}/${SAMPLE_ID}.transcript_models.gtf \
    ${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_ISO_asm.bed || exit 1

#create the file list which is the input for the merge step
printf "%s\tcapped\t1,1,1\tIsoQuant\n%s\tcapped\t1,1,1\tStringTie\n" \
    "${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_ISO_asm.bed" \
    "${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_STR_asm.bed" \
    > ${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_asm_file_list.txt || exit 1

# Merge parameters

file_list="${MERGED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_asm_file_list.txt"
prefix="${SAMPLE_ID}_${FLOWCELL}_brain_asm_merged"
fivethresh=300
splicethresh=10
threethresh=300

${TAMA_BASE_DIR}/tama_merge.py \
# merge the two bed12 files
python -m \
tama.tama_merge \
    -f ${file_list} \
    -p ${prefix} \
    -a ${fivethresh} \
    -m ${splicethresh} \
    -z ${threethresh} \
    -d merge_dup

# move the output files of the tama merge tool into to the merged file list directory
mv ./${prefix}* ${MERGED_PAIRWISE_DIR}/ || exit 1

# convert the merged bed12 file back to gtf
python \
${TAMA_CONVERTER}/tama_convert_bed_gtf_ensembl_no_cds.py \
    ${MERGED_PAIRWISE_DIR}/${prefix}.bed \
    ${MERGED_PAIRWISE_DIR}/${prefix}.gtf || exit 1

# Final step use IsoQuant to quantify the merged assembly
isoquant.py \
    -t 60 \
    --reference ${REF_FASTA} \
    --transcript_quantification all \
    --gene_quantification all \
    --no_model_construction \
    --data_type assembly \
    --count_exons   \
    --bam ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
    --genedb ${MERGED_PAIRWISE_DIR}/${prefix}.gtf \
    -o ${ISOQUANT_OUT} || exit 1

#making the directory for stringtie output
mkdir -p ${MERGED_PAIRWISE_DIR}/gffcompare/

#making the annotated gtf and comparison with the ref transcriptome using gffcompare
gffcompare \
    -r /data/CARDPB/resources/hg38/gencode.v43.annotation.gtf \
    -o ${MERGED_PAIRWISE_DIR}/gffcompare/${SAMPLE_ID} \
    ${MERGED_PAIRWISE_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_asm_merged.gtf

#Quantifying the annoted gtf from gffcompare with IsoQuant
isoquant.py \
    -t 60 \
    --reference ${REF_FASTA} \
    --transcript_quantification all \
    --gene_quantification all \
    --no_model_construction \
    --data_type assembly \
    --count_exons   \
    --bam ${MAPPED_DIR}/${SAMPLE_ID}_${FLOWCELL}_brain_mapped.sorted.bam \
    --genedb ${MERGED_PAIRWISE_DIR}/gffcompare/${SAMPLE_ID}.annotated.gtf \
    --prefix ${SAMPLE_ID} \
    -o ${ISOQUANT_OUT} || exit 1
