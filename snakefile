#!/usr/bin/env python3
import pandas as pd
import os

wildcard_constraints:
    sample_id = "[A-Za-z0-9_\-]+",
    flowcell_id = "[A-Za-z0-9_]+"

class DotDict(dict):
    """DotDict class allows accessing dictionary keys as attributes."""
    def __getattr__(self, attr):
        if attr in self:
            return self[attr]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{attr}'")
    def __setattr__(self, key, value):
        self[key] = value
    def __delattr__(self, item):
        try:
            del self[item]
        except KeyError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{item}'")


configfile: 'config.yml'
config = DotDict(config)

samples = pd.read_csv(config.sample_file, header=None, sep=None, engine='python')
#samples = pd.read_csv('basecalling_sample_sheet_jan_19.txt', header=None, sep=None, engine='python')

sample_flowcell_pairs = [{"sample_id": row[0], "flowcell_id": row[1]} for row in samples.values]

all_sample_names = samples.iloc[:,0].tolist()
all_flowcell_ids = samples.iloc[:,1].tolist()

# this `all` rule defines the final outputs at the very end of the workflow that need to be produced.
# Snakemake will then start thinking backwards to determine which rules are necessary
# to generate those final outputs.


rule all:
    input:
        expand(
            config.merge_dir + '/{sample_id}/{sample_id}_{flowcell_id}.annotated.gtf',
            zip,  # ensure pairing of sample_id and flowcell_id from each line
            sample_id   = samples.iloc[:,0].tolist(),
            flowcell_id = samples.iloc[:,1].tolist()
        )

rule basecall:
    input:
        pod5 = config.base_dir + '/{sample_id}/{sample_id}/{flowcell_id}/pod5'
    output:
        ubam = config.ont_ubam + '/{sample_id}/{sample_id}_{flowcell_id}.bam'
    resources:
        runtime=7200, mem_mb=150000, gpu=2, gpu_model='a100', disk_mb=50000
    threads: 30
    params:
        model = config.dorado_model,
        outdir = config.ont_ubam + '/{sample_id}/'
    shell:
        """
        
        scripts/basecalling.sh \
            --infile {input.pod5} \
            --outfile {output.ubam} \
            --outdir {params.outdir} \
            --model {params.model} \
        """

rule trimming:
    input:  
        ubam = config.ont_ubam + '/{sample_id}/{sample_id}_{flowcell_id}.bam'
    output:  
        fastq = config.pychopper_dir + '/{sample_id}/{sample_id}_{flowcell_id}.trimmed.fastq',
    params:
        outdir = config.pychopper_dir
    threads: 60
    resources:
        runtime=4320, mem_mb=120000, disk_mb=50000
    shell: 
        """
        scripts/trimming.sh \
            --infile {input.ubam} \
            --outfile {output.fastq}  \
            --outdir {params.outdir}
        """

rule alignment:
    input:  
        fastq = config.pychopper_dir + '/{sample_id}/{sample_id}_{flowcell_id}.trimmed.fastq'
    output:  
        mapped_bam = config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain_mapped.sorted.bam'
    threads: 120
    resources:
        runtime=4320, mem_mb=200000, disk_mb=100000
    params:
        mapped_dir = config.mapping_dir + '/{sample_id}/',
        human_fasta = config.genome,
        sirv_fasta = config.sirvome
    shell: 
        """
        scripts/alignment.sh \
            --infile {input.fastq} \
            --outfile {output.mapped_bam} \
            --outdir {params.mapped_dir} \
            --genome {params.human_fasta} \
            --sirvome {params.sirv_fasta}
        """

rule stringtie:
    input:
        mapped_bam = config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain_mapped.sorted.bam'
    output:
        sirv_stringtie_gtf = config.stringtie_dir + '/{sample_id}/sirv/{sample_id}_{flowcell_id}_sirv.stringtie.gtf',
        brain_stringtie_gtf = config.stringtie_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain.stringtie.gtf'
    threads: 60
    params:
        human_ref_gtf = config.human_genedb,
        sirv_ref_gtf = config.sirv_genedb,
        sir_stringtie_dir = config.stringtie_dir + '/{sample_id}/sirv/',
        brain_stringtie_dir = config.stringtie_dir + '/{sample_id}/'
    resources:
        runtime=4320, mem_mb=120000, disk_mb=50000, slurm_partition='norm'
    shell: 
        """
        scripts/sirv_stringtie_assembly.sh \
            --infile {input.mapped_bam} \
            --outfile {output.sirv_stringtie_gtf} \
            --outdir {params.sir_stringtie_dir} \
            --ref_gtf {params.sirv_ref_gtf}
        scripts/brain_stringtie_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.brain_stringtie_dir} \
            --outfile {output.brain_stringtie_gtf} \
            --ref_gtf {params.human_ref_gtf} 
        """


rule isoquant:
    input:  
        mapped_bam = config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain_mapped.sorted.bam'
    output:
        brain_isoquant_gtf = config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.transcript_models.gtf',
    threads: 60
    params:
        prefix = '{sample_id}_{flowcell_id}',
        sirv_prefix = '{sample_id}_{flowcell_id}_sirv',
        human_fasta = config.genome,
        sirv_fasta = config.sirvome,
        human_ref_gtf = config.human_genedb,
        sirv_ref_gtf = config.sirv_genedb,
        sirv_isoquant_dir = config.isoquant_dir + '/{sample_id}/sirv/',
        brain_isoquant_dir = config.isoquant_dir + '/{sample_id}/'
    resources:
        runtime=4320, mem_mb=120000, disk_mb=50000
    shell: 
        """
        scripts/sirv_isoquant_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.sirv_isoquant_dir} \
            --sirvome {params.sirv_fasta} \
            --genedb {params.sirv_ref_gtf} \
            --prefix {params.sirv_prefix}
        scripts/brain_isoquant_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.brain_isoquant_dir} \
            --genome {params.human_fasta} \
            --genedb {params.human_ref_gtf} \
            --prefix {params.prefix} 
        """

rule merge:
    input:
        brain_isoquant_gtf = config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.transcript_models.gtf',
        brain_stringtie_gtf = config.stringtie_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain.stringtie.gtf',
        mapped_bam = config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_brain_mapped.sorted.bam'
    threads: 120
    params:
        prefix = '{sample_id}_{flowcell_id}',
        human_fasta = config.genome,
        human_ref_gtf = config.human_genedb,
        tama_scripts = config.tama_scripts,
        merge_dir = config.merge_dir + '/{sample_id}',
    resources:
        runtime=4320, mem_mb=400000, disk_mb=100000
    output:
        annotated_gtf = config.merge_dir + '/{sample_id}/{sample_id}_{flowcell_id}.annotated.gtf'

    shell: 
        """
        scripts/transcript_merge.sh \
            --input_bam {input.mapped_bam} \
            --isoquant_gtf {input.brain_isoquant_gtf} \
            --stringtie_gtf {input.brain_stringtie_gtf} \
            --merge_dir {params.merge_dir} \
            --genome {params.human_fasta} \
            --ref_gtf {params.human_ref_gtf} \
            --prefix {params.prefix} \
            --outfile {output.annotated_gtf} \
            --tama_script_dir {params.tama_scripts} \
        """

