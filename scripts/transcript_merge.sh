#!/usr/bin/env bash

usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

            # Your options go here.
            --input_bam) assert_argument "$1" "$opt"; INPUT_BAM="$1"; shift;;
            --isoquant_gtf) assert_argument "$1" "$opt"; ISOQUANT_GTF="$1"; shift;;
            --stringtie_gtf) assert_argument "$1" "$opt"; STRINGTIE_GTF="$1"; shift;;
            --outfile) assert_argument "$1" "$opt"; OUTFILE="$1"; shift;;
            --genome) assert_argument "$1" "$opt"; GENOME="$1"; shift;;
            --merge_dir) assert_argument "$1" "$opt"; MERGE_DIR="$1"; shift;;
            --tama_script_dir) assert_argument "$1" "$opt"; TAMA_BASE_DIR="$1"; shift;;
            --ref_gtf) assert_argument "$1" "$opt"; REF_GTF="$1"; shift;;
            --prefix) assert_argument "$1" "$opt"; PREFIX="$1"; shift;;
            --threads) assert_argument "$1" "$opt"; THREADS="$1"; shift;;
            --isoquant_counts) assert_argument "$1" "$opt"; ISOQUANT_COUNTS="$1"; shift;;
            --assembly_mode) assert_argument "$1" "$opt"; ASSEMBLY_MODE="$1"; shift;;
            --scripts_dir) assert_argument "$1" "$opt"; SCRIPTS_DIR="$1"; shift;;
      
            # Arguments processing. You may remove any unneeded line after the 1st.
            -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
            -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
            *)          usage_error "this should NEVER happen ($opt)";;               # sanity test for previous patterns
    
        esac
    done
    shift  # $EOL
fi



# Rest of code

#make the output directory if it does not exist

mkdir -p ${MERGE_DIR}

#preprocess stringtie gtf to remove the . from strand column
#this is necessary for the merge step to work properly
awk -F '\t' '$7 != "." {print}' ${STRINGTIE_GTF} \
    > ${STRINGTIE_GTF%.gtf}_tmp.gtf

#convert stringtie gtf to bed12 (we don't need this file)
# python \
# ${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_stringtie.py \
#     ${STRINGTIE_GTF%.gtf}_tmp.gtf \
#     ${OUTFILE%.annotated.gtf}_stringtie.bed

# # # after you convert the stringtie gtf to bed12, you can remove the tmp gtf
# rm ${STRINGTIE_GTF%.gtf}_tmp.gtf || exit 1


# #convert isoquant to bed12 (we also don't need this file)
# python \
# ${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_ensembl.py \
#     ${ISOQUANT_GTF} \
#     ${OUTFILE%.annotated.gtf}_isoquant.bed


###############################################################################
# FILTERING BLOCK (discovery mode only)
# In discovery mode, filter both arms before TAMA merge:
#   1. Filter StringTie GTF by TPM > 0
#   2. Re-quantify filtered StringTie transcripts with IsoQuant (#2)
#   3. Filter IsoQuant direct arm by CPM > 1
#   4. Filter IsoQuant-on-StringTie arm by CPM > 1
#   5. Re-convert filtered GTFs to BED for TAMA merge
###############################################################################

#### We don't need these variables anymore since we will be using the filtered versions
# STRINGTIE_BED="${OUTFILE%.annotated.gtf}_stringtie.bed"
# ISOQUANT_BED="${OUTFILE%.annotated.gtf}_isoquant.bed"

if [ "${ASSEMBLY_MODE}" = "discovery" ]; then
    echo "=== Discovery mode: applying filtering before TAMA merge ==="

    # --- Step 1: Filter StringTie GTF to keep only transcripts with TPM > 0 ---
    echo "Step 1: Filtering StringTie GTF by TPM > 0..."
    STRINGTIE_FILTERED_GTF="${MERGE_DIR}/${PREFIX}_stringtie_tpm_filtered.gtf"
    python ${SCRIPTS_DIR}/filter_stringtie_by_tpm.py \
        --gtf ${STRINGTIE_GTF%.gtf}_tmp.gtf \
        --output-gtf ${STRINGTIE_FILTERED_GTF} \
        --tpm-threshold 0 || exit 1

    # # after you filter the stringtie processed gtf, you can remove the tmp gtf
    rm ${STRINGTIE_GTF%.gtf}_tmp.gtf || exit 1

    # --- Step 2: IsoQuant #2 — re-quantify filtered StringTie transcripts ---
    echo "Step 2: Running IsoQuant re-quantification on filtered StringTie GTF..."
    ISOQUANT_ON_STRINGTIE_DIR="${MERGE_DIR}/${PREFIX}_isoquant_on_stringtie"
    isoquant.py \
        -t "${THREADS}" \
        --reference ${GENOME} \
        --transcript_quantification unique_only \
        --gene_quantification unique_splicing_consistent \
        --no_model_construction \
        --data_type nanopore \
        --count_exons \
        --bam ${INPUT_BAM} \
        --genedb ${STRINGTIE_FILTERED_GTF} \
        --prefix ${PREFIX}_stringtie_requant \
        -o ${ISOQUANT_ON_STRINGTIE_DIR} || exit 1

    # Locate the IsoQuant-on-StringTie outputs
    ISOQUANT_ON_ST_COUNTS="${ISOQUANT_ON_STRINGTIE_DIR}/${PREFIX}_stringtie_requant/${PREFIX}_stringtie_requant.transcript_counts.tsv"
    
    
    ####GTF will not be produced for requantification so reuse the STRINGTIE_FILTERED_GTF
    # ISOQUANT_ON_ST_GTF="${ISOQUANT_ON_STRINGTIE_DIR}/${PREFIX}_stringtie_requant/${PREFIX}_stringtie_requant.transcript_models.gtf"


    # --- Step 3: Filter IsoQuant direct arm (#1) by CPM > 1 ---
    echo "Step 3: Filtering IsoQuant direct arm by CPM > 1..."
    ISOQUANT_CPM_FILTERED_GTF="${MERGE_DIR}/${PREFIX}_isoquant_cpm_filtered.gtf"
    python ${SCRIPTS_DIR}/filter_transcripts_by_cpm.py \
        --counts ${ISOQUANT_COUNTS} \
        --gtf ${ISOQUANT_GTF} \
        --output-gtf ${ISOQUANT_CPM_FILTERED_GTF} \
        --output-ids ${MERGE_DIR}/${PREFIX}_isoquant_passing_ids.txt \
        --cpm-threshold 1.0 || exit 1

    # --- Step 4: Filter IsoQuant-on-StringTie arm (#2) by CPM > 1 ---
    echo "Step 4: Filtering IsoQuant-on-StringTie arm by CPM > 1..."
    ####Change Isoquant_onST_gtf
    STRINGTIE_CPM_FILTERED_GTF="${MERGE_DIR}/${PREFIX}_stringtie_cpm_filtered.gtf"
    python ${SCRIPTS_DIR}/filter_transcripts_by_cpm.py \
        --counts ${ISOQUANT_ON_ST_COUNTS} \
        --gtf ${STRINGTIE_FILTERED_GTF} \
        --output-gtf ${STRINGTIE_CPM_FILTERED_GTF} \
        --output-ids ${MERGE_DIR}/${PREFIX}_stringtie_passing_ids.txt \
        --cpm-threshold 1.0 || exit 1

    # --- Step 5: Re-convert filtered GTFs to BED12 for TAMA merge ---
    echo "Step 5: Converting filtered GTFs to BED12..."

    # StringTie arm (was re-quantified by IsoQuant, so use ensembl converter)
    python ${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_ensembl.py \
        ${STRINGTIE_CPM_FILTERED_GTF} \
        ${OUTFILE%.annotated.gtf}_stringtie_filtered.bed || exit 1

    # IsoQuant direct arm
    python ${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_ensembl.py \
        ${ISOQUANT_CPM_FILTERED_GTF} \
        ${OUTFILE%.annotated.gtf}_isoquant_filtered.bed || exit 1

    # Use filtered BEDs for TAMA merge
    STRINGTIE_BED="${OUTFILE%.annotated.gtf}_stringtie_filtered.bed"
    ISOQUANT_BED="${OUTFILE%.annotated.gtf}_isoquant_filtered.bed"

    echo "=== Filtering complete. Proceeding to TAMA merge with filtered transcripts ==="
fi

#create the file list which is the input for the merge step
printf "%s\tcapped\t1,1,1\tStringTie\n%s\tcapped\t1,1,1\tIsoQuant\n" \
    "${STRINGTIE_BED}" \
    "${ISOQUANT_BED}" \
    > ${MERGE_DIR}/${PREFIX}_transcript_merge_list.txt || exit 1

# Merge parameters

file_list="${MERGE_DIR}/${PREFIX}_transcript_merge_list.txt"
name="${PREFIX}_merged"
fivethresh=300
splicethresh=10
threethresh=300

# merge the two bed12 files
python \
${TAMA_BASE_DIR}/tama_merge2.py \
    -f ${file_list} \
    -p ${name} \
    -a ${fivethresh} \
    -m ${splicethresh} \
    -z ${threethresh} \
    -d merge_dup

# move the output files of the tama merge script from working directory into to the merged file list directory
mv ./${name}* ${MERGE_DIR}/

# convert the merged bed12 file back to gtf
python \
${TAMA_BASE_DIR}/tama_convert_bed_gtf_ensembl_no_cds.py \
    ${MERGE_DIR}/${name}.bed \
    ${MERGE_DIR}/${name}.gtf || exit 1

#making the comparison with gffcompare
gffcompare \
    -r ${REF_GTF} \
    -o ${MERGE_DIR}/${PREFIX} \
    ${MERGE_DIR}/${name}.gtf

# Quantifying the annoted gtf from gffcompare with IsoQuant
isoquant.py \
    -t "${THREADS}" \
    --reference ${GENOME} \
    --transcript_quantification unique_only \
    --gene_quantification unique_splicing_consistent \
    --no_model_construction \
    --data_type nanopore \
    --count_exons   \
    --bam ${INPUT_BAM} \
    --genedb ${OUTFILE} \
    --prefix ${PREFIX} \
    -o ${MERGE_DIR} || exit 1
