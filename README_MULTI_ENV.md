# ONT Long-read RNA Sequencing Pipeline (Multi-Environment)

## Overview

This is a **refactored, portable Snakemake-based bioinformatics pipeline** for processing Oxford Nanopore Technologies (ONT) long-read RNA sequencing data. The pipeline works across **three environments**:

- **Biowulf HPC** — Full pipeline from basecalling onward
- **Generic HPC (SLURM)** — Pipeline from trimming onward (pre-basecalled inputs)
- **Local Workstation** — Pipeline from trimming onward (pre-basecalled inputs)

**Key improvement**: A single codebase that auto-detects your environment and configures itself appropriately.

---

## Quick Start by Environment

### Biowulf Users
```bash
# 1. Clone the repo and navigate to snakemake_pipeline
cd snakemake_pipeline

# 2. Edit input_example.txt with your sample sheet (tab-separated: sample_id flowcell_id)
# 3. Optionally customize config/biowulf.yaml if paths differ

# 4. Run (auto-detects Biowulf)
./snakemake.sh

# Or explicitly specify:
./snakemake.sh biowulf
```

### Generic HPC Users
```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate lrrna

# 2. Customize config/hpc_slurm.yaml with your cluster's reference paths
# 3. Place pre-basecalled BAM files in: snakemake/ONT_UBAM/{sample_id}/
# 4. Edit input_example.txt with your sample sheet

# 5. Run (auto-detects SLURM)
./snakemake.sh

# Or explicitly specify:
./snakemake.sh slurm
```

### Local Workstation Users
```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate lrrna

# 2. Customize config/default.yaml with your reference paths
# 3. Place pre-basecalled BAM files in: snakemake/ONT_UBAM/{sample_id}/
# 4. Edit input_example.txt with your sample sheet

# 5. Run (auto-detects local)
./snakemake.sh

# Or explicitly specify:
./snakemake.sh default
```

---

## Architecture & Key Concepts

### Configuration Hierarchy

The pipeline uses a **three-level configuration system**:

1. **`config/default.yaml`** — Universal baseline with placeholder paths
   - Relative output directories (works out-of-the-box)
   - Placeholder reference paths (users must customize)
   - `skip_basecall: true` (default: start from trimming)

2. **`config/biowulf.yaml`** — Biowulf-specific overrides
   - Absolute paths to `/data/CARDPB/resources/`
   - `skip_basecall: false` (enable basecalling)
   - Override only what differs from default

3. **`config/hpc_slurm.yaml`** — Generic HPC template
   - Commented placeholders for users to fill in
   - Inherits `skip_basecall: true` from default

**How it works**: Snakemake loads `default.yaml`, then layers environment-specific overrides on top.

### Environment Auto-Detection

`snakemake.sh` automatically detects your environment:

```bash
# Check 1: Is $SLURM_CLUSTER_NAME set to "biowulf"?
→ Use biowulf profile & modules

# Check 2: Is sinfo command available (generic SLURM)?
→ Use slurm profile & check for snakemake/singularity

# Check 3: Neither found?
→ Use default (local) profile & require conda environment
```

Override detection manually:
```bash
./snakemake.sh biowulf    # Force Biowulf
./snakemake.sh slurm      # Force generic HPC
./snakemake.sh default    # Force local
```

### The `skip_basecall` Flag

**Biowulf (`skip_basecall: false`)**:
- Basecalling rule is **enabled**
- Pipeline expects POD5 files in: `{sample_id}/{sample_id}/{flowcell_id}/pod5/`
- Snakemake runs: `basecall → trimming → alignment → assembly → merge`

**Local/HPC (`skip_basecall: true`)**:
- Basecalling rule is **disabled**
- Pipeline expects pre-basecalled BAM files in: `snakemake/ONT_UBAM/{sample_id}/*.bam`
- Snakemake runs: `trimming → alignment → assembly → merge`

This allows the same pipeline code to work in different environments without modification.

---

## Detailed Setup Instructions

### 1. Clone & Navigate

```bash
git clone https://github.com/NIH-CARD/CARDlongread_ONT_long_read_RNA.git
cd CARDlongread_ONT_long_read_RNA/snakemake_pipeline
```

### 2. Install Dependencies

**Biowulf users** (uses modules):
```bash
# Just run snakemake.sh — it loads modules automatically
./snakemake.sh biowulf --dry-run
```

**Local/HPC users** (uses conda):
```bash
conda env create -f environment.yml
conda activate lrrna
```

### 3. Configure Reference Paths

#### Biowulf
Edit `config/biowulf.yaml` — paths are already set for `/data/CARDPB/resources/`

#### Generic HPC
Edit `config/hpc_slurm.yaml` and uncomment/customize:
```yaml
# Example for a generic HPC cluster
genome: '/shared/genomes/GRCh38.fa'
human_genedb: '/shared/references/gencode.v43.annotation.gtf'
sirvome: '/shared/spike_ins/SIRV_multi-fasta.fasta'
sirv_genedb: '/shared/spike_ins/SIRV.gtf'
sirv_ref_gtf: '/shared/spike_ins/SIRV_ref.gtf'
```

Ask your HPC administrator for the location of these files.

#### Local Workstation
Edit `config/default.yaml` with absolute paths:
```yaml
genome: '/home/user/references/GRCh38.fa'
human_genedb: '/home/user/references/gencode.v43.gtf'
sirvome: '/home/user/references/SIRV.fasta'
sirv_genedb: '/home/user/references/SIRV_genedb.gtf'
sirv_ref_gtf: '/home/user/references/SIRV_ref.gtf'
```

### 4. Prepare Sample Sheet

Edit `input_example.txt` (tab-separated, no headers):
```
NABEC_KEN-1069_FTX_RNA	20250114_2319_2A_PAW72725_5119b730
NABEC_KEN-1092_FTX_RNA	20250114_2320_2B_PBA15382_a372bdf5
NABEC_SH-05-16_FTX_RNA	20250114_2315_3D_PBA14601_3935a955
```

Columns:
- **sample_id**: Unique identifier (e.g., `NABEC_KEN-1069_FTX_RNA`)
- **flowcell_id**: Sequencing run ID (e.g., `20250114_2319_2A_PAW72725_5119b730`)

### 5. Prepare Input Data

**Biowulf users**: POD5 files must be in:
```
{sample_id}/{sample_id}/{flowcell_id}/pod5/
```

**Local/HPC users**: Pre-basecalled BAM files must be in:
```
snakemake/ONT_UBAM/{sample_id}/{sample_id}_{flowcell_id}.bam
```

### 6. Dry-run (Validate Pipeline)

```bash
./snakemake.sh --dry-run
```

This shows what Snakemake will do without running anything. Useful for catching config errors early.

### 7. Run the Pipeline

```bash
# Local/HPC: runs in foreground
./snakemake.sh

# Biowulf: submit to cluster
sbatch snakemake.sh

# Or force a specific profile
./snakemake.sh biowulf --dry-run
```

---

## Pipeline Overview

### Rules & Workflow

```
POD5 files (Biowulf only)
    ↓
[basecall] — Convert POD5 → BAM using Dorado
    ↓
Pre-basecalled BAM or basecall output
    ↓
[trimming] — Trim/orient reads using Pychopper → FASTQ
    ↓
[alignment] — Map to GRCh38 + SIRV using Minimap2 → BAM
    ↓
[stringtie] — Assemble transcripts using StringTie → GTF
    ↓
[isoquant] — Assemble transcripts using IsoQuant → GTF
    ↓
[merge] — Merge StringTie + IsoQuant using TAMA → final GTF
    ↓
snakemake/ASSEMBLY/MERGED/{sample_id}/*.annotated.gtf
```

### Output Files

Final outputs per sample:
- `*.annotated.gtf` — Merged, reannotated transcript models
- Intermediate outputs in `snakemake/ASSEMBLY/`, `snakemake/MAPPED/`, etc.

### Resource Requirements

| Rule | CPU | Memory | GPU | Runtime |
|------|-----|--------|-----|---------|
| basecall | 30 | 150GB | 4×A100 | 72h |
| trimming | 120 | 120GB | — | 72h |
| alignment | 120 | 200GB | — | 72h |
| stringtie | 20 | 120GB | — | 72h |
| isoquant | 20 | 120GB | — | 72h |
| merge | 20 | 400GB | — | 72h |

---

## Configuration Reference

### Default Config Keys

```yaml
# Base directories
base_dir: './'                          # Root for relative paths
ont_ubam: 'snakemake/ONT_UBAM'         # Basecalled BAM files
pychopper_dir: 'snakemake/PYCHOPPER'   # Trimmed FASTQ
mapping_dir: 'snakemake/MAPPED'        # Aligned BAM
assembly_dir: 'snakemake/ASSEMBLY/'    # Assembly outputs
isoquant_dir: 'snakemake/ASSEMBLY/ISOQUANT'
stringtie_dir: 'snakemake/ASSEMBLY/STRINGTIE'
merge_dir: 'snakemake/ASSEMBLY/MERGED' # Final outputs

# Reference genomes
genome: '/path/to/GRCh38.fa'
human_genedb: '/path/to/gencode.v43.gtf'
sirvome: '/path/to/SIRV.fasta'
sirv_genedb: '/path/to/SIRV_genedb.gtf'
sirv_ref_gtf: '/path/to/SIRV_ref.gtf'

# Basecalling (Biowulf)
dorado_model: 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0'

# Control flow
skip_basecall: true|false              # Enable/disable basecalling

# Input
sample_file: 'input_example.txt'        # Sample sheet path
```

---

## Execution Examples

### Biowulf: Run Full Pipeline

```bash
# On login node
./snakemake.sh --dry-run

# When ready, submit to cluster
sbatch snakemake.sh
```

### Local: Dry-run & Debug

```bash
# Check what will run
./snakemake.sh --dry-run

# Run single rule (for testing)
snakemake trimming --profile ./snakemake_profiles/default

# View DAG
./snakemake.sh --dry-run --dag | dot -Tpdf > dag.pdf
```

### Generic HPC: SLURM Submission

```bash
# Interactive
srun ./snakemake.sh --dry-run
srun ./snakemake.sh

# Batch job
sbatch -N 1 -t 72:00:00 -J lrrna_pipeline snakemake.sh slurm
```

---

## Troubleshooting

### Basecalling Rule Not Running (Biowulf)

**Problem**: Pipeline starts at trimming instead of basecalling.

**Check**:
1. Verify `config/biowulf.yaml` has `skip_basecall: false`
2. Confirm POD5 files exist: `ls {sample_id}/{sample_id}/{flowcell_id}/pod5/`
3. Run: `./snakemake.sh --dry-run` to see rule DAG

### Pre-basecalled BAMs Not Found (Local/HPC)

**Problem**: Error: "Missing input file `snakemake/ONT_UBAM/...`"

**Check**:
1. BAM files must be in: `snakemake/ONT_UBAM/{sample_id}/{sample_id}_{flowcell_id}.bam`
2. Verify `input_example.txt` sample IDs match BAM filenames exactly
3. Ensure `skip_basecall: true` in your config

### Container Image Not Found

**Problem**: "Container image not found. Downloading..."

**Fix**: The script auto-downloads `lrrna_0.9.sif` on first run (requires internet). If offline, manually pull:
```bash
singularity pull oras://quay.io/datatecnica/lrrna:0.9
```

### Out of Memory on Merge Rule

**Problem**: Merge rule fails with OOM.

**Fix**: The merge rule requires 400GB memory. Check your HPC allocation:
```bash
# Biowulf
sbatch --mem=400G snakemake.sh

# Generic HPC
sbatch --mem=400G -t 72:00:00 snakemake.sh slurm
```

### Snakemake Not Found (Local/HPC)

**Problem**: "Error: snakemake not found."

**Fix**: Activate conda environment:
```bash
conda activate lrrna
./snakemake.sh
```

---

## Advanced Configuration

### Custom Thread Allocation

Edit `config/resources.yaml` to adjust CPU/memory per rule:
```yaml
resources:
  trimming:
    threads: 60          # Default 120
    memory_mb: 80000     # Default 120000
```

### Change Reference Version

To use GRCh39 instead of GRCh38:
1. Update `config/biowulf.yaml` (or your environment config):
   ```yaml
   genome: '/data/CARDPB/resources/hg39/GCA_000001405.26_GRCh39.fa'
   human_genedb: '/data/CARDPB/resources/hg39/gencode.v44.gtf'
   ```
2. Test with dry-run first
3. Monitor merge rule output for annotation compatibility

---

## Appendix: Data Structure Example (NABEC Cohort)

This section describes the internal data organization used for the NABEC cohort at NIH-CARD. It serves as an example of how to structure your data for the pipeline.

**Data Structure**    
Our raw ONT sequencing data are collected and transfered to their respective cohort directories. For example, our NABEC cohort has the following path: `/data/CARD_AUX/LRS_temp/NABEC_RNA` . Within the `NABEC_RNA` directory each sample folder are organized as follow:  

~~~  
/NABEC_RNA/NABEC_KEN-1069_FTX_RNA   
 └── NABEC_KEN-1069_FTX_RNA  
    └── 20250114_2319_2A_PAW72725_5119b730  
        ├── fastq_fail  
        ├── fastq_pass  
        ├── other_reports  
        └── pod5
~~~  
This file hierarchy allows us to easily retrieve the samples list used to run the pipeline scripts and `.json` file paths which we use to get sequencing summary statistics from [CARDlongread-report-parser](https://github.com/molleraj/CARDlongread-report-parser).  
  
**Generating Sample Sheets**  
Using the directory structure above we generate the samples list and `.json` file paths using the follow command and scripts.  

**1. Generate JSON report paths:**    
```bash
find /data/CARD_AUX/LRS_temp/NABEC_RNA/ \
-type f  \
-name "*.json" \
> /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/json_report_paths.txt  
```  

**2. Generate Sample Sheet (input_example.txt format):**
```bash
while read -r jsonfilepaths; do  
	SAMPLE_ID=$(basename "$(dirname "$(dirname "$jsonfilepaths")")")  
	FLOWCELL_ID=$(basename "$(dirname "$jsonfilepaths")")  
	echo -e "${SAMPLE_ID}\t${FLOWCELL_ID}"   \
   >>  /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/sample_sheet.txt
done < /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/json_report_paths.txt  
```  

**Output Example:**
```  
NABEC_KEN-1069_FTX_RNA  20250114_2319_2A_PAW72725_5119b730
NABEC_KEN-1092_FTX_RNA  20250114_2320_2B_PBA15382_a372bdf5
NABEC_KEN-1127_FTX_RNA  20250114_2320_2C_PAY66952_ea90759d
...
```

---

## File Structure

```
snakemake_pipeline/
├── README.md                         # Original documentation
├── README_MULTI_ENV.md              # This file
├── snakefile                         # Main workflow
├── snakemake.sh                      # Auto-detecting launcher
├── environment.yml                   # Conda dependencies
├── input_example.txt                 # Sample sheet (edit this)
│
├── config/
│   ├── default.yaml                 # Universal baseline
│   ├── biowulf.yaml                 # Biowulf overrides
│   ├── hpc_slurm.yaml               # Generic HPC template
│   └── resources.yaml                # Thread/memory specs
│
├── snakemake_profiles/
│   ├── default/config.yaml          # Local execution
│   ├── biowulf/config.yaml          # Biowulf SLURM profile
│   └── slurm/config.yaml            # Generic SLURM profile
│
├── scripts/
│   ├── basecalling.sh
│   ├── trimming.sh
│   ├── alignment.sh
│   ├── human_stringtie_assembly.sh
│   ├── sirv_stringtie_assembly.sh
│   ├── human_isoquant_assembly.sh
│   ├── sirv_isoquant_assembly.sh
│   ├── transcript_merge.sh
│   └── *.py                          # Helper scripts
│
└── logs/                             # Execution logs (auto-created)
```

---

## Contributing & Support

For issues or questions:
- Check the [original README.md](README.md) for detailed tool documentation
- Review config files for inline comments
- Test with `--dry-run` before submitting large jobs
- Contact: [NIH-CARD](https://github.com/NIH-CARD/CARDlongread_ONT_long_read_RNA)

---

## References

- **Dorado**: https://github.com/nanoporetech/dorado
- **Pychopper**: https://github.com/nanoporetech/pychopper
- **Minimap2**: https://github.com/lh3/minimap2
- **StringTie**: https://github.com/gpertea/stringtie
- **IsoQuant**: https://ablab.github.io/IsoQuant/
- **TAMA**: https://github.com/GenomeRIK/tama
- **gffcompare**: https://github.com/gpertea/gffcompare
