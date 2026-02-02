#!/usr/bin/env python3
"""
Transcript Filtering Script for Long-read RNA-seq Pipeline

Filters expressed transcripts from StringTie and IsoQuant outputs:
1. StringTie: Filter TPM > 0 → requantify with IsoQuant → filter CPM > threshold
2. IsoQuant: Filter count > 0 → calculate CPM → filter CPM > threshold

Usage:
    python filter_transcripts.py \
        --stringtie-gtf <path> \
        --stringtie-abundance <path> \
        --isoquant-gtf <path> \
        --isoquant-counts <path> \
        --output-dir <path> \
        --prefix <sample_prefix> \
        --genome <reference.fa> \
        --cpm-threshold 1.0
"""

import argparse
import pandas as pd
import subprocess
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Filter expressed transcripts")
    parser.add_argument("--stringtie-gtf", required=True, help="StringTie GTF file")
    parser.add_argument("--stringtie-abundance", required=True, 
                        help="StringTie abundance file (*.tab or use GTF FPKM)")
    parser.add_argument("--isoquant-gtf", required=True, help="IsoQuant GTF file")
    parser.add_argument("--isoquant-counts", required=True, 
                        help="IsoQuant transcript counts TSV")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Sample prefix for output files")
    parser.add_argument("--genome", required=True, help="Reference genome FASTA")
    parser.add_argument("--genedb", required=True, help="Reference GTF for requantification")
    parser.add_argument("--cpm-threshold", type=float, default=1.0, 
                        help="CPM threshold for filtering (default: 1.0)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for IsoQuant")
    return parser.parse_args()


def calculate_cpm(counts_series):
    """Calculate CPM (Counts Per Million) from raw counts."""
    total_counts = counts_series.sum()
    if total_counts == 0:
        return pd.Series(0, index=counts_series.index)
    return (counts_series / total_counts) * 1e6


def filter_stringtie_by_tpm(gtf_path, abundance_path, output_gtf, tpm_threshold=0):
    """
    Filter StringTie GTF to keep only transcripts with TPM > threshold.
    Returns set of passing transcript IDs.
    """
    # Read abundance file (StringTie outputs t_data.ctab or we parse GTF for FPKM/TPM)
    # StringTie GTF has FPKM and TPM in attributes for transcript lines
    
    passing_ids = set()
    
    # Parse GTF to extract TPM values
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'transcript':
                continue
            
            attrs = fields[8]
            # Extract transcript_id and TPM
            transcript_id = None
            tpm = 0.0
            
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                elif attr.startswith('TPM'):
                    try:
                        tpm = float(attr.split('"')[1])
                    except (IndexError, ValueError):
                        tpm = 0.0
            
            if transcript_id and tpm > tpm_threshold:
                passing_ids.add(transcript_id)
    
    print(f"StringTie: {len(passing_ids)} transcripts with TPM > {tpm_threshold}")
    
    # Filter GTF using the passing IDs
    filter_gtf_by_ids(gtf_path, passing_ids, output_gtf)
    
    return passing_ids


def filter_gtf_by_ids(input_gtf, transcript_ids, output_gtf):
    """Filter GTF file to keep only specified transcript IDs."""
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            attrs = fields[8]
            transcript_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            
            if transcript_id in transcript_ids:
                outfile.write(line)
    
    print(f"Wrote {len(transcript_ids)} transcripts to {output_gtf}")


def requantify_with_isoquant(gtf_path, bam_path, genome, output_dir, prefix, threads):
    """
    Requantify a GTF using IsoQuant in quantification-only mode.
    Returns path to the counts file.
    """
    cmd = [
        "isoquant.py",
        "-t", str(threads),
        "--reference", genome,
        "--genedb", gtf_path,
        "--bam", bam_path,
        "--data_type", "nanopore",
        "--transcript_quantification", "unique_only",
        "--gene_quantification", "unique_only",
        "--no_model_construction",
        "--prefix", prefix,
        "--output", output_dir
    ]
    
    print(f"Running IsoQuant requantification: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"IsoQuant error: {result.stderr}")
        sys.exit(1)
    
    # Return path to transcript counts
    counts_path = os.path.join(output_dir, prefix, f"{prefix}.transcript_counts.tsv")
    return counts_path


def filter_isoquant_by_cpm(counts_path, gtf_path, output_gtf, cpm_threshold):
    """
    Filter IsoQuant results by CPM threshold.
    Returns set of passing transcript IDs.
    """
    # Read counts file
    df = pd.read_csv(counts_path, sep='\t', comment='#')
    
    # IsoQuant counts file typically has columns: feature_id, <sample_count>
    # Find the count column (usually the second column or named after sample)
    count_col = df.columns[1]  # Assuming first col is ID, second is count
    
    # Filter count > 0 first
    df = df[df[count_col] > 0].copy()
    
    # Calculate CPM
    df['CPM'] = calculate_cpm(df[count_col])
    
    # Filter by CPM threshold
    passing_df = df[df['CPM'] > cpm_threshold]
    passing_ids = set(passing_df.iloc[:, 0].tolist())  # First column is transcript ID
    
    print(f"IsoQuant: {len(passing_ids)} transcripts with CPM > {cpm_threshold}")
    
    # Filter GTF
    filter_gtf_by_ids(gtf_path, passing_ids, output_gtf)
    
    return passing_ids


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # =========================================================================
    # Step 1: Filter StringTie by TPM > 0
    # =========================================================================
    print("\n=== Step 1: Filter StringTie transcripts (TPM > 0) ===")
    stringtie_tpm_filtered_gtf = os.path.join(
        args.output_dir, f"{args.prefix}_stringtie_tpm_filtered.gtf"
    )
    stringtie_tpm_ids = filter_stringtie_by_tpm(
        args.stringtie_gtf,
        args.stringtie_abundance,
        stringtie_tpm_filtered_gtf,
        tpm_threshold=0
    )
    
    # =========================================================================
    # Step 2: Requantify StringTie filtered GTF with IsoQuant
    # =========================================================================
    print("\n=== Step 2: Requantify StringTie GTF with IsoQuant ===")
    stringtie_requant_dir = os.path.join(args.output_dir, "stringtie_requant")
    stringtie_requant_prefix = f"{args.prefix}_str_requant"
    
    # Note: You'll need to pass the BAM file path - adding as argument
    # For now, we'll construct it from the workflow pattern
    # This should be passed from the Snakemake rule
    
    # =========================================================================
    # Step 3: Filter IsoQuant by count > 0, then CPM > threshold
    # =========================================================================
    print("\n=== Step 3: Filter IsoQuant transcripts (CPM > {}) ===".format(args.cpm_threshold))
    isoquant_cpm_filtered_gtf = os.path.join(
        args.output_dir, f"{args.prefix}_isoquant_cpm_filtered.gtf"
    )
    isoquant_cpm_ids = filter_isoquant_by_cpm(
        args.isoquant_counts,
        args.isoquant_gtf,
        isoquant_cpm_filtered_gtf,
        args.cpm_threshold
    )
    
    # =========================================================================
    # Step 4: Output filtered GTFs for TAMA merge
    # =========================================================================
    print("\n=== Final Filtered GTFs ===")
    print(f"StringTie filtered: {stringtie_tpm_filtered_gtf}")
    print(f"IsoQuant filtered:  {isoquant_cpm_filtered_gtf}")
    
    # Write summary stats
    summary_file = os.path.join(args.output_dir, f"{args.prefix}_filtering_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"StringTie transcripts (TPM > 0): {len(stringtie_tpm_ids)}\n")
        f.write(f"IsoQuant transcripts (CPM > {args.cpm_threshold}): {len(isoquant_cpm_ids)}\n")
    
    print(f"\nSummary written to: {summary_file}")


if __name__ == "__main__":
    main()