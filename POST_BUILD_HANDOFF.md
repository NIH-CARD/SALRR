# `lrrna_1.2.sif` — Post-Build Handoff

Container was built in a GitHub Codespace and is ready to ship.

| Item | Value |
|---|---|
| Image file | `lrrna_1.2.sif` (repo root) |
| Size | ~1.6 GB |
| Base | `docker://ubuntu:24.04` (glibc 2.39) |
| Env | conda env `lrrna` at `/opt/mamba/envs/lrrna/` (auto-activates) |
| Recipe | `singularity.def` (committed in `5ee3ef6`) |
| Branch | `feature/isomatch_replace_tama` |

**Smoke test passed as non-root**: bedtools 2.31.1, rustqc, isomatch 0.6.1, samtools 1.24, IsoQuant 4.0.0, stringtie 2.1.7, minimap2 2.31, gffcompare 0.12.10, multiqc 1.35, mosdepth 0.3.14, NanoPlot 1.46.2, pychopper — all work with plain `singularity exec lrrna_1.2.sif <tool>`.

---

## 1. Get the SIF off the Codespace

The image is 1.6 GB — do **NOT** `git add` it (already covered by `.gitignore` for 1.1; add 1.2 too — see step 5).

### Option A — GitHub Release asset (recommended, keeps a public copy)
```bash
# From the Codespace, on the repo root:
gh release create lrrna-container-1.2 \
  --repo NIH-CARD/CARDlongread_ONT_long_read_RNA \
  --title "lrrna container v1.2" \
  --notes "Adds bedtools, rustqc, IsoMatch v0.6.1. Base bumped to ubuntu:24.04 (glibc 2.39). Renames isoquant.py -> isoquant." \
  lrrna_1.2.sif
```
Anyone (incl. Cory / Biowulf) can then grab it with:
```bash
gh release download lrrna-container-1.2 \
  --repo NIH-CARD/CARDlongread_ONT_long_read_RNA \
  --pattern 'lrrna_1.2.sif'
```

### Option B — Direct download in Codespace terminal
```bash
# In Codespace:
python3 -m http.server 8000    # then use the forwarded port URL to download in browser
```
Save file locally, then `scp` to Biowulf.

### Option C — scp straight from Codespace (if you have Biowulf SSH set up here)
```bash
scp lrrna_1.2.sif <user>@biowulf.nih.gov:/data/CARD/<path>/
```

---

## 2. Pull to Biowulf

```bash
# On Biowulf, in the project data dir (NOT $HOME — SIFs are big):
cd /data/CARD/<your-project-path>
gh release download lrrna-container-1.2 \
  --repo NIH-CARD/CARDlongread_ONT_long_read_RNA \
  --pattern 'lrrna_1.2.sif'
ls -lh lrrna_1.2.sif

# Sanity check on a compute node (interactive):
sinteractive --mem=4g
module load singularity   # or apptainer, depending on Biowulf module
singularity exec lrrna_1.2.sif bash -lc 'bedtools --version && isomatch --version && isoquant --version'
```
Expect `Activating lrrna` + all three version lines.

---

## 3. Hand off to Cory for `quay.io/datatecnica/lrrna:1.2`

Cory has push rights to `quay.io/datatecnica/lrrna`. He needs to run (one-time login + push):
```bash
# On any machine with singularity + the SIF:
singularity remote login --username <quay-robot-user> oras://quay.io
singularity push lrrna_1.2.sif oras://quay.io/datatecnica/lrrna:1.2
```
Verify from anywhere:
```bash
singularity pull oras://quay.io/datatecnica/lrrna:1.2
```

---

## 4. Update pipeline scripts (do AFTER the image is on quay.io)

The new container renames the IsoQuant entry point: `isoquant.py` → `isoquant`. Six call sites to update:

- [scripts/human_isoquant_assembly.sh](scripts/human_isoquant_assembly.sh#L49)
- [scripts/sirv_isoquant_assembly.sh](scripts/sirv_isoquant_assembly.sh#L43)
- [scripts/transcript_merge.sh](scripts/transcript_merge.sh#L101)
- [scripts/transcript_merge.sh](scripts/transcript_merge.sh#L203)
- [scripts/transcript_merge_isomatch.sh](scripts/transcript_merge_isomatch.sh#L81)
- [scripts/transcript_merge_isomatch.sh](scripts/transcript_merge_isomatch.sh#L178)

One-shot sed (dry-run first):
```bash
# Dry-run — shows what would change:
grep -rn 'isoquant\.py' scripts/
# Apply:
sed -i 's/\bisoquant\.py\b/isoquant/g' \
  scripts/human_isoquant_assembly.sh \
  scripts/sirv_isoquant_assembly.sh \
  scripts/transcript_merge.sh \
  scripts/transcript_merge_isomatch.sh
# Verify:
grep -rn 'isoquant\.py' scripts/   # should print nothing
grep -rn '\bisoquant \\' scripts/  # sanity: should show the 6 lines
```

### Bump the image tag in the pipeline

**snakemake.sh line 66** — change `lrrna:1.1` to `lrrna:1.2`:
```bash
sed -i 's|oras://quay.io/datatecnica/lrrna:1.1|oras://quay.io/datatecnica/lrrna:1.2|' snakemake.sh
```

**snakefile** — 7 rules reference `./lrrna_1.1.sif`:
```bash
sed -i 's|./lrrna_1.1.sif|./lrrna_1.2.sif|g' snakefile
grep -n 'lrrna_1\.' snakefile   # verify all are now 1.2
```

---

## 5. `.gitignore`

Add the new SIF so it never accidentally gets committed:
```bash
grep -q '^lrrna_1.2.sif$' .gitignore || echo 'lrrna_1.2.sif' >> .gitignore
```

---

## 6. Quick end-to-end test on Biowulf

```bash
cd /data/CARD/<pipeline-workdir>
# Ensure lrrna_1.2.sif is present (from step 2) OR let snakemake.sh pull it via oras://
bash snakemake.sh --dry-run          # verify DAG builds
bash snakemake.sh --until human_isoquant_assembly   # one-rule smoke test on the test_sample/
```
Watch for:
- `Activating lrrna` in the log (env is auto-activating)
- No `isoquant.py: command not found` errors
- IsoMatch step (in `transcript_merge_isomatch.sh`) reaches its `isomatch classify` call

---

## 7. Commit + PR

```bash
git add scripts/human_isoquant_assembly.sh scripts/sirv_isoquant_assembly.sh \
        scripts/transcript_merge.sh scripts/transcript_merge_isomatch.sh \
        snakefile snakemake.sh .gitignore POST_BUILD_HANDOFF.md
git commit -m "Bump container to lrrna:1.2 and rename isoquant.py -> isoquant"
git push
gh pr create --fill --base main
```

---

## Appendix — known non-issues

- `WARNING: passwd file doesn't exist in container, not updating` on Codespace exec — Codespace-only cosmetic; does not appear on Biowulf.
- `WARNING: integrity: signature not found for object group 1` during build — SIF is unsigned; expected.
- `NanoPlot` prints a Kaleido/Chrome warning before its version — harmless; only affects static plot export.

## Appendix — recipe notes for future maintainers

The Codespace build hit a permission bug where `/.singularity.d/env/90-environment.sh` came out mode 754 (root-only readable), so non-root users' shells couldn't source the conda activation. Worked around by sandbox+chmod+repack. Permanent fix for next build:
```bash
sudo singularity build --fix-perms lrrna_1.2.sif singularity.def
```
