# Snakemake README


Build a graph of the rule hierarchy:

```bash
./snakemake.sh --dry-run --rulegraph | dot -Tpdf > rules.pdf
```

```bash
./snakemake.sh --dry-run --dag | dot -Tpdf > dag.pdf
```