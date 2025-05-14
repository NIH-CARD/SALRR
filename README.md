# NIA CARD Long Read RNA Sequencing Pipeline
![Workflow](data_processing_workflow.png)

**Overview**

This repository provides a modular and comprehensive pipeline designed for processing Oxford Nanopore Technologies (ONT) long-read RNA sequencing data, optimized for human brain samples. It includes essential steps such as base calling, read trimming, rescue and reorientation, mapping, transcript assembly, and the merging of assembled transcripts.

Currently, the pipeline is implemented as a series of individual scripts, each handling a distinct step. The ultimate goal is to integrate these components into a robust Snakemake workflow to improve reproducibility, scalability, and ease of use; supporting the transcriptomics efforts of NIA CARD and the broader research community.

**Requirements**

Ensure the following dependencies and tools are installed before running the pipeline:
              *Dorado (base calling)
              *Pychopper (read trimming, rescue and orientation)
              *Minimap2 (Mapping)
              *IsoQuant(transcript assembly)
              *StringTie (transcript assembly)
              *TAMA (assembly merging)
Addition dependencies and tools:
              *CARD_Long_read_report_parser(sequencing report parsing and visualization)
              *SIRVsuite(SIRV RNA Spike-in Control QC)

Instructions for installing each tool are provided in their respective documentation:

- **[Dorado_Basecalling](https://github.com/nanoporetech/dorado)**
- **[Pychopper](https://github.com/nanoporetech/pychopper)**
- **[Minimap2](https://github.com/lh3/minimap2)**
- **[Stringtie](https://github.com/gpertea/stringtie) and [IsoQuant](https://ablab.github.io/IsoQuant/index.html)**
- **[TAMA](https://github.com/GenomeRIK/tama.git)**

Additional tools:
- **[CARD_Long_read_report_parser](https://github.com/molleraj/CARDlongread-report-parser)** 
- **[SIRVsuite](https://github.com/Lexogen-Tools/SIRVsuite)**

**Pipeline Steps**
1. Base Calling
We perform base calling on ONT `.pod5` files using Dorado v9 to generate `.bam` files.
Example Command:
`sbatch --array=1-4 dorado_v090_basecalling_array_RNA.sh`

