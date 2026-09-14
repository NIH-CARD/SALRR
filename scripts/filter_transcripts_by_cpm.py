#!/usr/bin/env python3
"""
Filter transcripts by CPM threshold from IsoQuant count files.

Usage:
    python filter_transcripts_by_cpm.py \
        --counts <isoquant_transcript_counts.tsv> \
        --gtf <input.gtf> \
        --output-gtf <filtered.gtf> \
        --output-ids <passing_ids.txt> \
        --cpm-threshold 1.0
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter transcripts by CPM from IsoQuant counts"
    )
    parser.add_argument(
        "--counts", required=True,
        help="IsoQuant transcript counts TSV file"
    )
    parser.add_argument(
        "--gtf", required=True,
        help="Input GTF to filter"
    )
    parser.add_argument(
        "--output-gtf", required=True,
        help="Output filtered GTF"
    )
    parser.add_argument(
        "--output-ids", required=True,
        help="Output file with passing transcript IDs"
    )
    parser.add_argument(
        "--cpm-threshold", type=float, default=1.0,
        help="CPM threshold for filtering (default: 1.0)"
    )
    return parser.parse_args()


def read_counts(counts_path):
    """Read IsoQuant transcript counts, return dict of {transcript_id: count}."""
    counts = {}
    with open(counts_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            transcript_id = fields[0]
            try:
                count = float(fields[1])
            except ValueError:
                continue
            counts[transcript_id] = count
    return counts


def filter_by_cpm(counts, cpm_threshold):
    """Calculate CPM and return set of IDs passing threshold."""
    total = sum(counts.values())
    if total == 0:
        print("WARNING: Total counts = 0, no transcripts pass filter",
              file=sys.stderr)
        return set()

    passing = set()
    for tid, count in counts.items():
        cpm = (count / total) * 1e6
        if cpm > cpm_threshold:
            passing.add(tid)

    print(f"Total assigned reads: {total:.0f}", file=sys.stderr)
    print(f"CPM threshold: {cpm_threshold}", file=sys.stderr)
    print(f"Transcripts passing: {len(passing)} / {len(counts)} "
          f"(count > 0: {sum(1 for c in counts.values() if c > 0)})",
          file=sys.stderr)
    return passing


def filter_gtf(input_gtf, passing_ids, output_gtf):
    """Filter GTF keeping only lines with transcript_id in passing_ids."""
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

            # Extract transcript_id from attributes
            attrs = fields[8]
            transcript_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    # Handle both transcript_id "X" and transcript_id X
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        transcript_id = parts[1]
                    else:
                        transcript_id = attr.split()[-1]
                    break

            if transcript_id and transcript_id in passing_ids:
                outfile.write(line)
                kept += 1

    print(f"GTF lines: {kept} / {total} kept", file=sys.stderr)


def main():
    args = parse_args()

    # Read counts
    counts = read_counts(args.counts)

    # Filter by CPM
    passing_ids = filter_by_cpm(counts, args.cpm_threshold)

    # Write passing IDs
    with open(args.output_ids, 'w') as f:
        for tid in sorted(passing_ids):
            f.write(f"{tid}\n")

    # Filter GTF
    filter_gtf(args.gtf, passing_ids, args.output_gtf)


if __name__ == "__main__":
    main()