#!/usr/bin/env bash

usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
require_cmd () { command -v "$1" >/dev/null 2>&1 || usage_error "required command not found in PATH: $1"; }
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
            --isomatch_bin) assert_argument "$1" "$opt"; ISOMATCH_BIN="$1"; shift;;
      
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

mkdir -p "${MERGE_DIR}"

require_cmd bedtools

###############################################################################
# FILTERING BLOCK (discovery mode only)
# In discovery mode, filter both arms before IsoMatch merge:
#   1. Filter StringTie GTF by TPM > 0
#   2. Re-quantify filtered StringTie transcripts with IsoQuant (#2)
#   3. Filter IsoQuant direct arm by CPM > 1
#   4. Filter IsoQuant-on-StringTie arm by CPM > 1
#   5. Sort + index filtered GTFs for IsoMatch merge
###############################################################################


# Default behavior (quantification): merge direct StringTie and IsoQuant outputs.
STRINGTIE_MERGE_INPUT="${STRINGTIE_GTF}"
ISOQUANT_MERGE_INPUT="${ISOQUANT_GTF}"

if [ "${ASSEMBLY_MODE}" = "discovery" ]; then
    echo "=== Discovery mode: applying filtering before IsoMatch merge ==="

    # --- Step 1: Filter StringTie GTF to keep only transcripts with TPM > 0 ---
    echo "Step 1: Filtering StringTie GTF by TPM > 0..."
    STRINGTIE_FILTERED_GTF="${MERGE_DIR}/${PREFIX}_stringtie_tpm_filtered.gtf"
    python ${SCRIPTS_DIR}/filter_stringtie_by_tpm.py \
        --gtf "${STRINGTIE_GTF}" \
        --output-gtf "${STRINGTIE_FILTERED_GTF}" \
        --tpm-threshold 0 || exit 1

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

    # Discovery uses filtered transcripts for merge.
    STRINGTIE_MERGE_INPUT="${STRINGTIE_CPM_FILTERED_GTF}"
    ISOQUANT_MERGE_INPUT="${ISOQUANT_CPM_FILTERED_GTF}"
else
    echo "=== Quantification mode: merging StringTie and IsoQuant outputs directly ==="
fi

# --- Step 5: Sort + index merge inputs for IsoMatch ---
echo "Step 5: Sorting and indexing GTFs for IsoMatch..."

STRINGTIE_SORTED_GTF="${MERGE_DIR}/${PREFIX}_stringtie_for_merge.sorted.gtf"
ISOQUANT_SORTED_GTF="${MERGE_DIR}/${PREFIX}_isoquant_for_merge.sorted.gtf"

bedtools sort -i ${STRINGTIE_MERGE_INPUT} > ${STRINGTIE_SORTED_GTF} || exit 1
bedtools sort -i ${ISOQUANT_MERGE_INPUT} > ${ISOQUANT_SORTED_GTF} || exit 1

# Resolve IsoMatch from explicit CLI/env override, local checkout, then PATH.
if [ -n "${ISOMATCH_BIN:-}" ]; then
    ISOMATCH="${ISOMATCH_BIN}"
elif [ -n "${ISOMATCH:-}" ]; then
    ISOMATCH="${ISOMATCH}"
elif [ -x "${SCRIPTS_DIR}/../isomatch/target/release/isomatch" ]; then
    ISOMATCH="${SCRIPTS_DIR}/../isomatch/target/release/isomatch"
else
    ISOMATCH="isomatch"
fi

if ! command -v "${ISOMATCH}" >/dev/null 2>&1; then
    usage_error "IsoMatch executable not found. Set --isomatch_bin /path/to/isomatch or export ISOMATCH=/path/to/isomatch"
fi

echo "Using IsoMatch executable: ${ISOMATCH}"
"${ISOMATCH}" --version || usage_error "failed to execute IsoMatch binary: ${ISOMATCH}"

"${ISOMATCH}" index --ref-fa ${GENOME} ${STRINGTIE_SORTED_GTF} || exit 1
"${ISOMATCH}" index --ref-fa ${GENOME} ${ISOQUANT_SORTED_GTF} || exit 1

echo "=== Proceeding to IsoMatch merge ==="

# --- Step 6: Merge GTFs with IsoMatch ---
MERGE_PREFIX="${MERGE_DIR}/${PREFIX}_merged"
"${ISOMATCH}" merge --ref-fa ${GENOME} -o ${MERGE_PREFIX} ${STRINGTIE_SORTED_GTF} ${ISOQUANT_SORTED_GTF} || exit 1

# --- Step 7: Run IsoMatch classify on merged output ---
"${ISOMATCH}" index --ref-fa ${GENOME} ${MERGE_PREFIX}.merged.gtf.gz || exit 1
"${ISOMATCH}" classify \
    --ref-gtf ${REF_GTF} \
    --ref-fa ${GENOME} \
    --out ${MERGE_PREFIX} \
    ${MERGE_PREFIX}.merged.gtf.gz || exit 1

# Decompress classify output to the Snakemake-declared final output path.
zcat ${MERGE_PREFIX}.annotated.gtf.gz > ${OUTFILE} || exit 1

# --- Step 8: Re-quantify annotated merged GTF with IsoQuant ---
isoquant.py \
    -t "${THREADS}" \
    --reference ${GENOME} \
    --transcript_quantification unique_only \
    --gene_quantification unique_splicing_consistent \
    --no_model_construction \
    --data_type nanopore \
    --count_exons \
    --bam ${INPUT_BAM} \
    --genedb ${OUTFILE} \
    --prefix ${PREFIX} \
    -o ${MERGE_DIR} || exit 1
