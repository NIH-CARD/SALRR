# Snakemake README

#NOTE: First download snakemake profile
#if on biowulf, for example:

```bash
if [[ ! -d snakemake_profile ]]; then
    git clone https://github.com/NIH-HPC/snakemake_profile.git
fi
```

Build a graph of the rule hierarchy:

```bash
./snakemake.sh --dry-run --rulegraph | dot -Tpdf > rules.pdf
```

```bash
./snakemake.sh --dry-run --dag | dot -Tpdf > dag.pdf
```