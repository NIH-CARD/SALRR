#!/usr/bin/env python3
"""Aggregate per-sample IsoQuant outputs into cohort-wide matrices."""

import argparse
import glob
import os
import sys
import pandas as pd

# build a matrix from isoquant output file types

FILE_TYPES = {"transcript_counts": "transcript_counts_matrix.tsv",
    "transcript_tpm":    "transcript_tpm_matrix.tsv",
    "gene_counts":       "gene_counts_matrix.tsv",
    "gene_tpm":          "gene_tpm_matrix.tsv",
}

# create argument parsers

def parse_args():
    p = argparse.ArgumentParser(
        description="Build cohort count/TPM matrices from per-sample IsoQuant outputs."
    )
    p.add_argument("--quant-dir", required=True,
                   help="Directory containing per-sample IsoQuant output subfolders.")
    p.add_argument("--matrix-dir", required=True,
                   help="Output directory for the cohort matrices.")
    return p.parse_args()

# function to build the matrix

def build_matrix(quant_dir, file_type):
    """Join every sample's <file_type> file into one matrix (rows=feature_id, cols=sample)."""
    paths = sorted(glob.glob(os.path.join(quant_dir, "*", f"*.{file_type}.tsv")))
    if not paths:
        sys.exit(f"ERROR: no *.{file_type}.tsv files found under {quant_dir}")

    per_sample = []
    for path in paths:
        sample = os.path.basename(path).replace(f".{file_type}.tsv", "")
        col = pd.read_csv(path, sep="\t")
        col.columns = ["feature_id", sample]
        per_sample.append(col.set_index("feature_id")[sample])

    matrix = pd.concat(per_sample, axis=1).fillna(0)
    matrix.index.name = "feature_id"
    return matrix

def main():
    args = parse_args()
    os.makedirs(args.matrix_dir, exist_ok=True)
    for file_type, out_name in FILE_TYPES.items():
        matrix = build_matrix(args.quant_dir, file_type)
        out_path = os.path.join(args.matrix_dir, out_name)
        matrix.to_csv(out_path, sep="\t")
        print(f"Wrote {out_path}  ({matrix.shape[0]} features x {matrix.shape[1]} samples)")

if __name__ == "__main__":
    main()