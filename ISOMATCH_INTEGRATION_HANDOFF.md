# IsoMatch Integration Handoff (Snakemake Pipeline)

## Purpose
This document captures what is now implemented for IsoMatch-based per-sample merge in the Snakemake pipeline, and what remains for the next phase.

## Repository and Branch
- Repo path: `/data/CARD_AUX/LRS_temp/NABEC_RNA/snakemake_pipeline_isomatch`
- Active branch: `feature/isomatch_replace_tama`

## High-Level Goal
Integrate IsoMatch into the per-sample merge stage and use IsoMatch classify outputs as the annotated transcript model input for final re-quantification.

## Design Decisions Made
1. Use `isomatch classify` annotation output as the final per-sample annotated GTF (`--outfile`) instead of using `gffcompare`.
2. Keep discovery-specific filtering logic (TPM/CPM filtering) before merge.
3. In quantification mode, skip discovery filtering but still run merge, classify, and final IsoQuant re-quantification.
4. Keep IsoMatch binary resolution flexible: explicit CLI path, env var, local checkout binary, then PATH fallback.

## What Was Implemented

### 1) New merge script created
- New file: `scripts/transcript_merge_isomatch.sh`
- Script is tracked and committed.

### 2) Discovery filtering retained
In discovery mode, script performs:
1. StringTie TPM filtering
2. IsoQuant re-quant on filtered StringTie
3. IsoQuant-arm CPM filtering
4. StringTie-arm CPM filtering
5. Uses filtered GTFs as merge inputs

### 3) Quantification mode behavior aligned
- Quantification mode now merges directly from:
   - StringTie GTF
  - IsoQuant GTF
- No discovery filtering in quantification mode.

### 4) IsoMatch merge path implemented
Common steps (both modes):
1. Sort merge inputs with `bedtools sort`
2. Index each sorted input with `isomatch index`
3. Merge with `isomatch merge --ref-fa ${GENOME}`

### 5) IsoMatch classify added
After merge:
1. Index merged query GTF
2. Run `isomatch classify` against reference GTF and FASTA

### 6) Downstream contract now uses classify output
After classify:
1. Decompress `MERGE_PREFIX.annotated.gtf.gz` into `OUTFILE` (Snakemake-declared per-sample final annotated GTF)
2. Run final IsoQuant re-quantification with `--genedb ${OUTFILE}`

### 7) Snakefile wiring completed
- Merge rule now calls `scripts/transcript_merge_isomatch.sh`.

## Current Workflow State (Important)

### Current per-sample merge behavior
1. Merge rule runs per sample/flowcell pair.
2. For each sample pair, merge uses that sample's StringTie and IsoQuant outputs only.
3. Classify output is used to generate that sample's final annotated GTF (`*.annotated.gtf`).

## Verification Already Performed
1. Script-level syntax check passed:
   - `bash -n scripts/transcript_merge_isomatch.sh`
2. Merge rule updated to call `scripts/transcript_merge_isomatch.sh`.
3. Changes were committed and pushed on branch `feature/isomatch_replace_tama`.

## Remaining Work (Next Session Checklist)

### A) Run dry-run validation on updated wiring
1. Run:
   - `ml snakemake`
   - `snakemake -n --config profile=default --cores 1`
2. Confirm DAG resolves and includes expected per-sample IsoMatch merge path.

### B) Optional runtime validation (small scope)
1. Run one sample/subset test if possible.
2. Confirm expected key outputs exist:
   - merged GTF (`*.merged.gtf.gz`)
   - classify output (`*.classification.txt.gz`)
   - classify annotated GTF output (`*.annotated.gtf.gz`)
   - final IsoQuant requant outputs

### C) Improve install UX for IsoMatch (recommended)
1. Add optional config key (for example `isomatch_bin`) in `config/default.yaml`.
2. Pass that key into merge rule and forward to script via `--isomatch_bin`.
3. Document three supported modes in README:
   - explicit `--isomatch_bin`
   - `ISOMATCH` env var
   - PATH/local fallback

### D) New feature request: optional cross-sample merge and re-quantification
Add a second, optional workflow path where users can merge across many samples (for example ~40 samples) and then re-quantify against that cohort-level merged annotation.

Suggested design for this feature:
1. Keep existing per-sample merge rule as default behavior.
2. Add optional cohort mode toggle in config (for example `enable_cohort_merge: true/false`).
3. Add a new rule to collect per-sample annotated outputs (or per-sample merge inputs) into a cohort merge list.
4. Run `isomatch merge` across the cohort list, then `isomatch classify` on cohort merged GTF.
5. Produce cohort final annotated GTF (for example `snakemake/ASSEMBLY/MERGED/cohort/cohort.annotated.gtf`).
6. Add optional cohort re-quant rule to run IsoQuant with cohort annotated GTF for each sample BAM.
7. Allow running this cohort stage independently after an initial per-sample pipeline run.

Inputs for "run later" cohort mode:
1. Existing per-sample BAM files (already produced)
2. Existing per-sample StringTie/IsoQuant outputs or per-sample annotated outputs
3. Same reference FASTA + reference GTF

## Known Caveats / Notes
1. `isomatch classify` handling in prior testing indicated indexing and input format are important; merged query should be indexed before classify.
2. Tool availability depends on environment/modules (e.g., Snakemake loaded via `ml snakemake`).
3. Current script resolves IsoMatch binary as:
   - `--isomatch_bin` value (if provided)
   - `ISOMATCH` env var (if set)
   - `${SCRIPTS_DIR}/../isomatch/target/release/isomatch`
   - fallback to `isomatch` in `PATH`
4. Git pull does not install IsoMatch; users still need a runnable IsoMatch binary available by one of the supported resolution methods.

## Ready-to-Paste Prompt for Next AI Session
Use this if you want to quickly resume:

"Continue IsoMatch integration in `/data/CARD_AUX/LRS_temp/NABEC_RNA/snakemake_pipeline_isomatch` on branch `feature/isomatch_replace_tama`. Per-sample merge is now wired to `scripts/transcript_merge_isomatch.sh`; classify output is decompressed to `--outfile` and used for final IsoQuant requant. Next, run dry-run validation (`ml snakemake && snakemake -n --config profile=default --cores 1`) and then design optional cohort-level cross-sample merge + classify + re-quant that can run either inline or after per-sample outputs already exist." 