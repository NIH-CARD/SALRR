import argparse
import os
import polars as pl
from gtfparse import read_gtf

# Get input arguments from SLURM script
parser = argparse.ArgumentParser(description="Extract arguments from GTF files")
parser.add_argument('gtf_file', help="Path to GTF file to parse")
parser.add_argument('sample_id', help="Sample ID from the sample sheet")
parser.add_argument('flowcell_id', help="Flowcell ID from the sample sheet")
args = parser.parse_args()

# Read GTF
df = read_gtf(args.gtf_file)

# Filter for transcripts
features = df.filter(pl.col("feature") == "transcript")[["reference_id", "ref_gene_id", "FPKM"]]

# Rename columns
transcript_counts = features.rename({
    "reference_id": "tracking_ID",
    "ref_gene_id": "gene_ID",
    "FPKM": "FPKM_CHN"
})

# Define output directory and filename
base_dir = "/data/CARDPB/data/LRS_RNA/data/NGD"
output_dir = f"{base_dir}/SIRV_analysis/{args.sample_id}/"
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f"{args.sample_id}_{args.flowcell_id}_transcript_counts.tsv")

# Save as tab-delimited file
transcript_counts.write_csv(output_file, separator="\t")

print(f"Processed {args.gtf_file} → {output_file}")

# Define SIRVsuite directories
sirv_sheet_dir = f"{base_dir}/SIRV_analysis/"
alignment_path = os.path.join(
    f"{base_dir}/MAPPED/{args.sample_id}/",
    f"{args.sample_id}_{args.flowcell_id}_SIRV_mapped_filtered.sorted.bam"
)

# Create a sample sheet for the SIRVsuite
sirvsuite_sample_sheet = os.path.join(sirv_sheet_dir, "sirvsuite_sample_sheet.csv")
header = "sample_name;alignment_path;counting_path;read_orientation;counting_feature;replication_group\n"

# Check if the file exists to avoid adding multiple headers
file_exists = os.path.isfile(sirvsuite_sample_sheet)
   
with open(sirvsuite_sample_sheet, "a", encoding="utf-8") as f:
    if not file_exists:
        f.write(header)  # Write header if the file does not exist
    f.write(f"{args.sample_id};{alignment_path};{output_file};FWD;transcript;none\n")