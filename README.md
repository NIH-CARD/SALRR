# NIA CARD Long Read RNA Sequencing Pipeline
![Workflow](data_processing_workflow.png)

**Overview**

This repository provides a modular and comprehensive pipeline designed for processing Oxford Nanopore Technologies (ONT) long-read RNA sequencing data, optimized for human brain samples. It includes essential steps such as base calling, read trimming, rescue and reorientation, mapping, transcript assembly, and the merging of assembled transcripts.

Currently, the pipeline is implemented as a series of individual scripts, each handling a distinct step. The ultimate goal is to integrate these components into a robust Snakemake workflow to improve reproducibility, scalability, and ease of use; supporting the transcriptomics efforts of NIA CARD and the broader research community.

**Requirements**

Ensure the following dependencies and tools are installed before running the pipeline:  
   - Dorado (base calling)  
   - Pychopper (read trimming, rescue and orientation)  
   - Minimap2 (Mapping)  
   - IsoQuant(transcript assembly)  
   - StringTie (transcript assembly)  
   - TAMA , GFFCOMPARE, ISOQUANT (assembly merging, reannotation, and requantification)  

Addition dependencies and tools:  
   - CARD_Long_read_report_parser(sequencing report parsing and visualization)   
   - SIRVsuite(SIRV RNA Spike-in Control QC)  

Instructions for installing each tool are provided in their respective documentation:  

- **[Dorado_Basecalling](https://github.com/nanoporetech/dorado)**
- **[Pychopper](https://github.com/nanoporetech/pychopper)**
- **[Minimap2](https://github.com/lh3/minimap2)**
- **[Stringtie](https://github.com/gpertea/stringtie) and [IsoQuant](https://ablab.github.io/IsoQuant/index.html)**
- **[TAMA](https://github.com/GenomeRIK/tama.git) and [GFFCOMPARE](https://github.com/gpertea/gffcompare.git)**

Additional tools:
- **[CARD_Long_read_report_parser](https://github.com/molleraj/CARDlongread-report-parser)** 
- **[SIRVsuite](https://github.com/Lexogen-Tools/SIRVsuite)**

**Pipeline Steps**
1. **Base Calling**   
We perform base calling on ONT `.pod5` files using Dorado v9 to generate `.bam` files.  
Example Command:  
`sbatch --array=1-4 pipeline/01_basecalling/dorado_v090_basecalling_array_RNA.sh`  
2. **Read Trimming**  
Identify, orient, and trim full-length Nanopore cDNA reads using Pychopper. `.bam` files are converted to `.fastq` at the start of the run. 
Example Command:  
`sbatch --array=1-4 pipeline/02_trimming/pychopper.sh`  
3. **Mapping**  
Map SIRV RNA spike-in control reads to the SIRVome and extracting full-length brain sample reads to the reference genome GRCh38 using Minimap2  
Example Command:  
`sbatch --array=1-4 pipeline/03_alignment/minimap_brain_sirv.sh`  
4. **Assembly**  
Perform transcripts assembly using IsoQuant and StringTie on brain mapped reads and SIRV RNA spike-in control reads    
- Brain example command:  
`sbatch --array=1-4 pipeline/04_assembly/brain_assembly.sh`  
- SIRV example command:  
`sbatch --array=1-4 pipeline/04_assembly/sirv_assembly.sh`  
5. **Merging**  
Merging transcript assemblies from IsoQuant and StringTie using TAMA to generate a unified `.GTF` file, reannotation the `.GTF` with gffcompare and requantifying using IsoQuant
Example command:  
`sbatch --array=1-4 pipeline/05_merge/tama_merge.sh`  
  
**Output**  
The main final outputs of the pipeline are the merged `.GTF` file of the sample including the `.tsv` files with transcripts and gene quantifications.
# Data Structure and Pipeline Usage 
**Data Struture**    
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
This file hierarchy allows us to easily retrieve the samples list used to run the pipeline scripts and `.json` file paths which we use to get sequencing summary stastics from [CARDlongread-report-parser](https://github.com/molleraj/CARDlongread-report-parser).  
  
**Pipeline Usage**  
Using the directory structure above we extract make the samples list and `.json` file paths using thee follow command and scripts.  
- json file paths example command:    
```
find /data/CARD_AUX/LRS_temp/NABEC_RNA/ \  
-type f  \  
-name "*.json" \  
> /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/json_report_paths.txt  
```  
output:json_report_paths.txt    
```  
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_KEN-1069_FTX_RNA/NABEC_KEN-1069_FTX_RNA/20250114_2319_2A_PAW72725_5119b730/report_PAW72725_20250114_2321_5119b730.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_KEN-1092_FTX_RNA/NABEC_KEN-1092_FTX_RNA/20250114_2320_2B_PBA15382_a372bdf5/report_PBA15382_20250114_2326_a372bdf5.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_KEN-1127_FTX_RNA/NABEC_KEN-1127_FTX_RNA/20250114_2320_2C_PAY66952_ea90759d/report_PAY66952_20250114_2322_ea90759d.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_KEN-1142_FTX_RNA/NABEC_KEN-1142_FTX_RNA/20250114_2321_2D_PAY64818_eea1f032/report_PAY64818_20250114_2327_eea1f032.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-00-34_FTX_RNA/NABEC_SH-00-34_FTX_RNA/20250114_2321_2E_PAY64520_2d2a8c92/report_PAY64520_20250114_2327_2d2a8c92.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-03-17_FTX_RNA/NABEC_SH-03-17_FTX_RNA/20250116_1843_3G_PAY64563_27732d10/report_PAY64563_20250116_1850_27732d10.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-05-16_FTX_RNA/NABEC_SH-05-16_FTX_RNA/20250114_2315_3D_PBA14601_3935a955/report_PBA14601_20250114_2321_3935a955.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-97-09_FTX_RNA/NABEC_SH-97-09_FTX_RNA/20250114_2314_3B_PAW72691_179afafe/report_PAW72691_20250114_2320_179afafe.json
/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-97-53_FTX_RNA/NABEC_SH-97-53_FTX_RNA/20250114_2315_3C_PAY68370_975be5bd/report_PAY68370_20250114_2321_975be5bd.json  
```  
- samples list example command:  
```
while read -r jsonfilepaths; do  
	SAMPLE_ID=$(basename "$(dirname "$(dirname "$jsonfilepaths")")")  
	FLOWCELL_ID=$(basename "$(dirname "$jsonfilepaths")")  
	echo -e "${SAMPLE_ID}\t${FLOWCELL_ID}"   \
   >>  /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/sample_sheet.txt
done < /data/CARD_AUX/LRS_temp/NABEC_RNA/SAMPLE_SHEETS/json_report_paths.txt  
```  
output: sample_sheet.txt   
```  
NABEC_KEN-1069_FTX_RNA  20250114_2319_2A_PAW72725_5119b730
NABEC_KEN-1092_FTX_RNA  20250114_2320_2B_PBA15382_a372bdf5
NABEC_KEN-1127_FTX_RNA  20250114_2320_2C_PAY66952_ea90759d
NABEC_KEN-1142_FTX_RNA  20250114_2321_2D_PAY64818_eea1f032
NABEC_SH-00-34_FTX_RNA  20250114_2321_2E_PAY64520_2d2a8c92
NABEC_SH-03-17_FTX_RNA  20250116_1843_3G_PAY64563_27732d10
NABEC_SH-05-16_FTX_RNA  20250114_2315_3D_PBA14601_3935a955
NABEC_SH-97-09_FTX_RNA  20250114_2314_3B_PAW72691_179afafe
NABEC_SH-97-53_FTX_RNA  20250114_2315_3C_PAY68370_975be5bd  
```
With our samples list ready we can running the pipeline scripts sequentially for 4 samples.
```
# Base Calling  
sbatch --array=1-4 pipeline/01_basecalling/dorado_v090_basecalling_array_RNA.sh  

# Trimming  
sbatch --array=1-4 pipeline/02_trimming/pychopper.sh  

# Mapping    
sbatch --array=1-4 pipeline/03_alignment/minimap_brain_sirv.sh  

# Assembly   
sbatch --array=1-4 pipeline/04_assembly/brain_assembly.sh # mapping read to GRCh38   
sbatch --array=1-4 pipeline/04_assembly/sirv_assembly.sh  # mapping reads to SIRV genome  

# Merging    
sbatch --array=1-4 pipeline/05_merge/tama_merge_sh   
```  

