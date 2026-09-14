# NIA CARD Long-read RNA Sequencing Pipeline

![Workflow](data_processing_workflow.png)

A Snakemake pipeline for processing Oxford Nanopore Technologies (ONT) long-read
RNA sequencing data: read trimming, alignment, transcript assembly,
quantification, and optional transcript merging. It runs on the NIH Biowulf
cluster, generic SLURM clusters, and local workstations through a single
codebase with selectable execution profiles.

> **Note on basecalling:** Basecalling is currently supported on Biowulf only.
> On generic SLURM clusters and local workstations, the pipeline starts from
> pre-basecalled (unaligned) BAM files.

---

## Environments and Profiles

Select an execution profile when launching the pipeline:

```bash
./snakemake.sh <profile> [snakemake options]
```

| Profile   | Environment          | Starting point                     |
|-----------|----------------------|------------------------------------|
| `biowulf` | NIH Biowulf cluster  | POD5 (basecalling) or BAM          |
| `slurm`   | Generic SLURM HPC    | Pre-basecalled BAM                 |
| `default` | Local workstation    | Pre-basecalled BAM                 |

The profile is your choice; each one loads the appropriate modules or
environment for that system.

---

## Requirements

On the host system you only need:

- **Snakemake** (workflow engine)
- **Singularity/Apptainer** (containerized execution)

On Biowulf these are loaded automatically by the launcher. All pipeline tools
are bundled in a container that is pulled automatically on first run:

```
oras://quay.io/datatecnica/lrrna:1.2  ->  lrrna_1.2.sif
```

### Tools in the Pipeline

Processing:

- [Dorado](https://github.com/nanoporetech/dorado): basecalling (Biowulf only)
- [Pychopper](https://github.com/nanoporetech/pychopper): read trimming, rescue, and orientation
- [Minimap2](https://github.com/lh3/minimap2): alignment
- [StringTie](https://github.com/gpertea/stringtie): transcript assembly
- [IsoQuant](https://ablab.github.io/IsoQuant/index.html): assembly and quantification
- [IsoMatch](https://github.com/zhengxinchang/isomatch): assembly merging and reannotation

Quality control and reporting:

- [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/): contamination screening
- [samtools](https://www.htslib.org/) / [mosdepth](https://github.com/brentp/mosdepth): alignment and coverage stats
- [cramino (nanopack)](https://github.com/wdecoster/nanopack): long-read summary stats
- [MultiQC](https://multiqc.info/): aggregated QC report

### Additional tools

- [RSeQC](https://rseqc.sourceforge.net/) and [RustQC](https://github.com/seqeralabs/rustqc): additional RNA-seq QC checks
- [CARDlongread-report-parser](https://github.com/molleraj/CARDlongread-report-parser): NIH CARD available resource for parsing and summarizing sequencing run output

---

## Quick Start

```bash
# 1. Clone and enter the pipeline directory
git clone https://github.com/NIH-CARD/CARDlongread_ONT_long_read_RNA.git
cd CARDlongread_ONT_long_read_RNA

# 2. Set reference paths in the config for your profile
#    biowulf -> config/biowulf.yaml
#    slurm   -> config/slurm.yaml
#    default -> config/default.yaml

# 3. Edit the sample sheet (input_example.txt)

# 4. Validate the workflow (dry run)
./snakemake.sh <profile> -n

# 5. Run
./snakemake.sh <profile>
```

On Biowulf, submit as a batch job:

```bash
sbatch snakemake.sh biowulf
```

---

## Assembly Modes

Set `assembly_mode` in your config file:

- **`discovery`**: assembles novel transcripts with IsoQuant and StringTie, then
  merges them with IsoMatch into a unified, reannotated GTF. Final per-sample
  output is an annotated GTF. This is the most common mode.
- **`quantification`**: reference-only estimation with no novel transcript
  models. Final per-sample outputs are IsoQuant transcript/gene count tables and
  the StringTie GTF. The merge step is not run in this mode.

---

## Input Data

### Biowulf (POD5, basecalling enabled)

POD5 files are expected under:

```
{sample_id}/{sample_id}/{flowcell_id}/pod5/
```

For example:

```text
/NABEC_RNA/NABEC_KEN-1069_FTX_RNA
└── NABEC_KEN-1069_FTX_RNA
  └── 20250114_2319_2A_PAW72725_5119b730
    ├── fastq_fail
    ├── fastq_pass
    ├── other_reports
    └── pod5
```

The pipeline uses the sample ID, flowcell ID, and `pod5` directory at the
flowcell level to locate the input data.

### Generic HPC / Local (pre-basecalled BAM)

Unaligned BAM files are expected under:

```
snakemake_test/ONT_UBAM/{sample_id}/{sample_id}_{flowcell_id}.bam
```

### FastQ Screen reference databases

FastQ Screen requires the reference databases to be prepared locally and indexed before use. The pipeline expects database paths in the FastQ Screen config file to point to the prebuilt reference base names, and users should follow the official FastQ Screen documentation for creating those indexes. The example config file in this repository is a template and will need to be updated to match the local database paths on your system.

### Sample Sheet

Provide a tab-separated file with one sample per line and no header:

```
NABEC_KEN-1069_FTX_RNA	20250114_2319_2A_PAW72725_5119b730
NABEC_SH-05-16_FTX_RNA	20250114_2315_3D_PBA14601_3935a955
```

- **Column 1** — sample ID
- **Column 2** — flowcell ID

If you do not have a flowcell ID, any unique run label works in the second
column (for example an experiment name), as long as it matches your input file
names.

---

## Outputs and QC

The workflow creates a few main output areas. A typical layout looks like this:

```text
snakemake_test/
├── MAPPED/
│   └── {sample_id}/
│       └── {sample_id}_{flowcell_id}_human_mapped.sorted.bam
├── ASSEMBLY/
│   ├── ISOQUANT/
│   │   └── {sample_id}/
│   │       └── {sample_id}_{flowcell_id}/
│   │           ├── {sample_id}_{flowcell_id}.transcript_counts.tsv
│   │           └── {sample_id}_{flowcell_id}.gene_counts.tsv
│   ├── STRINGTIE/
│   │   └── {sample_id}/
│   │       └── {sample_id}_{flowcell_id}_human.stringtie.gtf
│   └── MERGED/
│       └── {sample_id}/
│           └── {sample_id}_{flowcell_id}.annotated.gtf
├── QC/
│   ├── multiqc_report.html
│   ├── mapping_qc/
│   ├── fastq_screen/
│   └── trimming_qc/
└── ASSEMBLY/COHORT/   # only when cohort merge is enabled
    ├── gtf/
    │   └── {cohort_name}.annotated.gtf
    ├── quant/
    ├── matrix/
    │   ├── transcript_counts_matrix.tsv
    │   └── gene_counts_matrix.tsv
    └── cohort_manifest.tsv
```

The main outputs to check after a successful run are:

- the mapped BAM file for each sample
- the sample-level annotated GTF in `discovery` mode
- the sample-level transcript count table and StringTie GTF in `quantification` mode
- the final MultiQC report in the QC directory
- the cohort merged GTF and matrices only when `use_cross_sample_merge` is enabled

A run is typically considered successful when the expected files above are present, Snakemake finishes without failed jobs, and the MultiQC HTML report is generated.
