#!/usr/bin/env python3
"""
Filter StringTie GTF to keep only transcripts with TPM > threshold.

Usage:
    python filter_stringtie_by_tpm.py \
        --gtf <stringtie.gtf> \
        --output-gtf <filtered.gtf> \
        --tpm-threshold 0
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter StringTie GTF by TPM threshold"
    )
    parser.add_argument(
        "--gtf", required=True,
        help="Input StringTie GTF file"
    )
    parser.add_argument(
        "--output-gtf", required=True,
        help="Output filtered GTF"
    )
    parser.add_argument(
        "--tpm-threshold", type=float, default=0.0,
        help="TPM threshold (default: 0, keeps transcripts with TPM > 0)"
    )
    return parser.parse_args()


def get_passing_transcript_ids(gtf_path, tpm_threshold):
    """Parse StringTie GTF, return set of transcript IDs with TPM > threshold."""
    passing = set()
    total = 0

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue

            total += 1
            attrs = fields[8]
            transcript_id = None
            tpm = 0.0

            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        transcript_id = parts[1]
                elif attr.startswith('TPM'):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        try:
                            tpm = float(parts[1])
                        except ValueError:
                            tpm = 0.0

            if transcript_id and tpm > tpm_threshold:
                passing.add(transcript_id)

    print(f"StringTie transcripts with TPM > {tpm_threshold}: "
          f"{len(passing)} / {total}", file=sys.stderr)
    return passing


def filter_gtf(input_gtf, passing_ids, output_gtf):
    """Filter GTF keeping only lines whose transcript_id is in passing_ids."""
    kept = 0
    total = 0

    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            total += 1

            attrs = fields[8]
            transcript_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        transcript_id = parts[1]
                    break

            if transcript_id and transcript_id in passing_ids:
                outfile.write(line)
                kept += 1

    print(f"GTF lines: {kept} / {total} kept", file=sys.stderr)


def main():
    args = parse_args()

    # Get passing IDs from transcript lines
    passing_ids = get_passing_transcript_ids(args.gtf, args.tpm_threshold)

    # Filter full GTF (transcript + exon lines)
    filter_gtf(args.gtf, passing_ids, args.output_gtf)


if __name__ == "__main__":
    main()