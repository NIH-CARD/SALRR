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

# loading the modules
# module load python/3.12
# module load isoquant/3.6.2
# module load gffcompare/0.12.6

#make the output directory if it does not exist

mkdir -p ${MERGE_DIR}

#preprocess stringtie gtf to remove the . from strand column
#this is necessary for the merge step to work properly
awk -F '\t' '$7 != "." {print}' ${STRINGTIE_GTF} \
    > ${STRINGTIE_GTF%.gtf}_tmp.gtf

#convert stringtie gtf to bed12
python \
${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_stringtie.py \
    ${STRINGTIE_GTF%.gtf}_tmp.gtf \
    ${OUTFILE%.annotated.gtf}_stringtie.bed

# # after you convert the stringtie gtf to bed12, you can remove the tmp gtf
rm ${STRINGTIE_GTF%.gtf}_tmp.gtf || exit 1


#convert isoquant to bed12
python \
${TAMA_BASE_DIR}/tama_format_gtf_to_bed12_ensembl.py \
    ${ISOQUANT_GTF} \
    ${OUTFILE%.annotated.gtf}_isoquant.bed

#create the file list which is the input for the merge step
printf "%s\tcapped\t1,1,1\tIsoQuant\n%s\tcapped\t1,1,1\tStringTie\n" \
    "${OUTFILE%.annotated.gtf}_stringtie.bed" \
    "${OUTFILE%.annotated.gtf}_isoquant.bed" \
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
    -t 20 \
    --reference ${GENOME} \
    --transcript_quantification all \
    --gene_quantification all \
    --no_model_construction \
    --data_type assembly \
    --count_exons   \
    --bam ${INPUT_BAM} \
    --genedb ${OUTFILE} \
    --prefix ${PREFIX} \
    -o ${MERGE_DIR} || exit 1
