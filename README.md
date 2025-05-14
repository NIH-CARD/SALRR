# **NIA CARD Long Read RNA Sequencing Pipeline
![Workflow](CARDlongread_ONT_long_read_RNA/data_processing_workflow.png)

**Overview**

This repository provides a modular and comprehensive pipeline designed for processing Oxford Nanopore Technologies (ONT) long-read RNA sequencing data, optimized for human brain samples. It includes essential steps such as base calling, read trimming, rescue and reorientation, mapping, transcript assembly, and the merging of assembled transcripts.

Currently, the pipeline is implemented as a series of individual scripts, each handling a distinct step. The ultimate goal is to integrate these components into a robust Snakemake workflow to improve reproducibility, scalability, and ease of use; supporting the transcriptomics efforts of NIA CARD and the broader research community.


- **[LRS bioinformatics tutorial](https://github.com/molleraj/lrs-bioinformatics-tutorial)**
- **[CARD_Long_read_report_parser](https://github.com/molleraj/CARDlongread-report-parser)**
- **[Dorado_Basecalling](https://github.com/nanoporetech/dorado)**
- **[Pychopper](https://github.com/nanoporetech/pychopper)**
- **[Minimap2](https://github.com/lh3/minimap2)**
- **[SIRVsuite](https://github.com/Lexogen-Tools/SIRVsuite)**
- **[Stringtie](https://github.com/gpertea/stringtie) and [IsoQuant](https://ablab.github.io/IsoQuant/index.html)**
- **[rnaseqtools_gtfmerge](https://github.com/Kingsford-Group/rnaseqtools)**
    