#!/usr/bin/env python3
import pandas as pd
import os

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
all_sample_names = samples.iloc[:,0].tolist()


# this `all` rule defines the final outputs at the very end of the workflow that need to be produced.
# Snakemake will then start thinking backwards to determine which rules are necessary
# to generate those final outputs.

rule all:
    input:
        expand("assembly/{sample_id}.merged.gtf", sample_id=all_sample_names)
       #bam = expand("basecalling/{sample_id}.bam", sample_id=all_sample_names)

rule basecall:
    input:
        pod5 = 'rawdata/{sample_id}.pod5'
    output:
        bam = 'basecalling/{sample_id}.bam',
        summary = 'basecalling/{sample_id}.summary.txt'
    resources:
        runtime=2880, mem_mb=300000, gpu=1, gpu_model='v100x'
    params:
        model = config.dorado_model
    envmodules:
        'dorado/0.9.0',
        'pod5/0.3.15'
    shell:
        """
        scripts/basecalling.sh \
            --infile {input.pod5} \
            --outfile {output.bam} \
            --model {params.model}
        """

rule trimming:
    input:  'basecalling/{sample_id}.bam'
    output:  'trimming/{sample_id}.trimmed.fastq'
    shell: 'scripts/trimming.sh'

rule alignment:
    input:  'trimming/{sample_id}.trimmed.fastq'
    output:  'alignment/{sample_id}.bam'
    shell: 'scripts/alignment.sh'

rule isoquant:
    input:  'alignment/{sample_id}.bam'
    output:  'assembly/{sample_id}.isoquant.gtf'
    shell: 'scripts/assembly.sh'

rule stringtie:
    input:  'alignment/{sample_id}.bam'
    output:  'assembly/{sample_id}.stringtie.gtf'
    shell: 'scripts/assembly.sh'

rule merge:
    input:  'assembly/{sample_id}.stringtie.gtf',
            'assembly/{sample_id}.isoquant.gtf'
    output:  'assembly/{sample_id}.merged.gtf'
    shell: 'scripts/merge.sh'

