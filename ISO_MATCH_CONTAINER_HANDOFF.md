# IsoMatch Container Integration Handoff

## Goal

Build a reproducible container image for the LRNA pipeline that includes:
1. `bedtools`
2. `rustqc`
3. `IsoMatch`

Then wire the pipeline to use that container cleanly, so users do not need a separate preinstalled environment.

Biowulf is currently unavailable, so all work in this session should be done with local repository inspection, container recipe edits, and non-Biowulf validation only.

## Current Understanding

1. The pipeline already expects a container image named `lrrna_1.1.sif`.
2. The pipeline launcher auto-pulls that image from Quay in [`snakemake.sh`](snakemake.sh).
3. The workflow itself already uses that container in [`snakefile`](snakefile).
4. This repo does not appear to contain the container build recipe itself.
5. The container build logic likely lives in another repo or branch, but this repo’s current `environment.yml` and container recipe content were shared manually by the user.

## Relevant Repo Files in This Workspace

1. [`environment.yml`](environment.yml)
2. [`snakemake.sh`](snakemake.sh)
3. [`snakefile`](snakefile)
4. [`scripts/transcript_merge_isomatch.sh`](scripts/transcript_merge_isomatch.sh)
5. [`config/default.yaml`](config/default.yaml)
6. [`README_MULTI_ENV.md`](README_MULTI_ENV.md)

## Assumptions to Keep in Mind

1. `bedtools` is a standard package and should be added in the normal package layer.
2. `rustqc` is likely best added via a pinned release binary or a pinned source build.
3. `IsoMatch` should be treated the same way as `rustqc`: pinned and installed during image build, not at runtime.
4. The final image should be versioned explicitly, not treated as a floating latest build.
5. The pipeline should not assume Biowulf is available during development or testing.

## Recommended Installation Strategy

### bedtools
Add it in the conda dependency file or the container’s conda install step, because it is a standard bioinformatics package.

### rustqc
Prefer one of these:
1. Use a pinned release binary if the project publishes one.
2. If no release binary exists, build from a pinned tag or commit during image build.
3. If it is a Rust crate with Cargo support, use Cargo only if the build is reproducible and the version is pinned.

### IsoMatch
Prefer one of these:
1. Use a pinned release binary if available.
2. If the project is Rust-based and meant to be built from source, build it from a pinned tag or commit in the container recipe.
3. Avoid unpinned git `main` builds and avoid runtime download/install inside Snakemake jobs.

## Implementation Plan

### Phase 1: Find the container build source of truth
1. Identify the repo or branch that contains the container recipe used to build `lrrna`.
2. Confirm whether the recipe is Singularity, Apptainer, Docker, or a build pipeline that generates the SIF.
3. Confirm how the current image version `lrrna_1.1` was produced.

### Phase 2: Decide the exact install method for each tool
1. Add `bedtools` to the package dependency layer.
2. Decide whether `rustqc` should be installed from a release artifact or built from source.
3. Decide whether `IsoMatch` should be installed from a release artifact or built from source.
4. Pin exact versions or exact commit hashes for `rustqc` and `IsoMatch`.
5. Record the chosen versions in the recipe comments or release notes.

### Phase 3: Update the container recipe
1. Install `bedtools` in the package layer.
2. Add `rustqc` installation steps in the image build.
3. Add `IsoMatch` installation steps in the image build.
4. Make sure the installed binaries are placed on `PATH` inside the container.
5. Add any needed build dependencies for Rust-based compilation.
6. Clean up build artifacts to keep the image smaller.

### Phase 4: Add validation inside the container build
1. Run tool version checks as part of image build or a post-build smoke test.
2. Confirm `bedtools` runs.
3. Confirm `rustqc` runs.
4. Confirm `IsoMatch` runs.
5. Confirm the binary paths are the ones expected inside the image.

### Phase 5: Rebuild and tag the container
1. Build a new image tag instead of overwriting the old one silently.
2. Example pattern: `lrrna_1.2` or another explicit versioned tag.
3. Keep the old image available until the new one is validated.
4. Record the image tag in the repo docs.

### Phase 6: Wire the pipeline to the updated image
1. Update [`snakemake.sh`](snakemake.sh) to pull the new image tag.
2. Update [`snakefile`](snakefile) so all rules use the new image path.
3. Verify that the merge rule still calls [`scripts/transcript_merge_isomatch.sh`](scripts/transcript_merge_isomatch.sh) correctly.
4. If needed, add config knobs for tool paths or image path overrides, but keep the default path simple.

### Phase 7: Improve the pipeline behavior for missing tools
1. If the container image is missing, the launcher should give a clear error or pull it explicitly.
2. If a tool is missing from the container, fail early with a helpful message.
3. Avoid runtime guessing or silent fallback logic for IsoMatch once it is containerized.

## Validation Plan Without Biowulf

Because Biowulf is down, validation should be local and container-centered.

1. Validate the container recipe syntax if possible.
2. Build the container locally or in the environment available to the agent.
3. Smoke test the image by checking tool versions.
4. Verify the pipeline launcher points to the correct image.
5. Do a non-Biowulf dry-run only if Snakemake is available locally.
6. If Snakemake is not available locally, do static validation of the launcher, recipe, and merge script wiring instead.

## What to Pay Attention To in the Merge Script

The merge script should not try to install tools dynamically at runtime if the image already contains them.

It should:
1. Assume the container provides `bedtools`, `rustqc`, and `IsoMatch`.
2. Use the tools directly from `PATH` inside the container.
3. Keep the current merge and classify workflow intact.
4. Fail with a clear error if a required binary is missing.

## Acceptance Criteria

The work is done when:
1. `bedtools` is included in the container image.
2. `rustqc` is included in the container image.
3. `IsoMatch` is included in the container image.
4. The image is versioned and reproducible.
5. The pipeline points to the new image.
6. The merge script runs against the container-provided binaries.
7. Tool version smoke tests succeed.
8. The workflow wiring remains intact.

## Suggested Order of Work

1. Locate the container recipe source.
2. Identify the current image build method.
3. Add `bedtools`, `rustqc`, and `IsoMatch` to the recipe.
4. Rebuild the image.
5. Smoke test the image.
6. Update the pipeline to reference the new image tag.
7. Verify merge script wiring.
8. Document the new supported setup in the README.

## Notes for the Next Agent

1. Prefer minimal, pinned, reproducible changes.
2. Do not redesign the workflow unless required.
3. Do not rely on Biowulf for validation in this session.
4. If a tool has both binary and source-install options, prefer the most reproducible one.
5. If the container recipe is in another repository, locate that repo first before editing anything in the pipeline repo.