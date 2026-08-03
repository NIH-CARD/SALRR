#!/usr/bin/env python3
import pandas as pd
import os

wildcard_constraints:
    sample_id = "[A-Za-z0-9_.\\-]+",  # sample_id can contain alphanumeric characters, underscores, periods, and/or hyphens
    flowcell_id = "[A-Za-z0-9_]+" #replace with project_id, or add project id? do we need to keep flowcells separate? perhaps for QC metrics

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

#loading environment-specific config based PROFILE environment variable
environment = config.get('profile', 'default')
configfile: f'config/{environment}.yaml'
config = DotDict(config)

# Read skip_basecall setting (default: True for local/HPC, can be overriden to False for Biowulf)
SKIP_BASECALL = config.get('skip_basecall', True)

#Read Enable or Disable SIRV spike-in analysis setting (default: set to True)
USE_SIRV = config.get('use_sirv', True)
if USE_SIRV and not config.get('sirvome'):
    raise ValueError("use_sirv enabled but sirvome path not set in config")
if USE_SIRV and not config.get('sirv_genedb'):
    raise ValueError("use_sirv enabled but sirv_genedb path not set in config")

# Read Enable or Disable FastQ Screen contamination screening (default: True)
USE_FASTQ_SCREEN = config.get('use_fastq_screen', True)
if USE_FASTQ_SCREEN and not config.get('fastq_screen_conf'):
    raise ValueError("use_fastq_screen enabled but fastq_screen_conf path not set in config")

# Read assembly mode setting (default: quantification for reference transcript abundance)
ASSEMBLY_MODE = config.get('assembly_mode', 'quantification')
if ASSEMBLY_MODE not in ['discovery', 'quantification']:
    raise ValueError(f"assembly_mode must be 'discovery' or 'quantification', got: {ASSEMBLY_MODE}")

samples = pd.read_csv(config.sample_file, header=None, sep=None, engine='python')

sample_flowcell_pairs = [{"sample_id": row[0], "flowcell_id": row[1]} for row in samples.values]

all_sample_names = samples.iloc[:,0].tolist()
all_flowcell_ids = samples.iloc[:,1].tolist()

# this `all` rule defines the final outputs at the very end of the workflow that need to be produced.
# Snakemake will then start thinking backwards to determine which rules are necessary
# to generate those final outputs.



rule all:
    input:
        expand(
            config.base_dir + config.merge_dir + '/{sample_id}/{sample_id}_{flowcell_id}.annotated.gtf',
            zip,  # ensure pairing of sample_id and flowcell_id from each line
            sample_id   = samples.iloc[:,0].tolist(),
            flowcell_id = samples.iloc[:,1].tolist()
        ),
        # Adding MultiQC report as a target to run after all samples are processed
        config.base_dir + config.qc_dir + '/multiqc_report.html'

# Basecalling rule only runs on Biowulf (skip_basecall: false)
if not SKIP_BASECALL:
    rule basecall:
        input:
            pod5 = config.base_dir + '{sample_id}/{sample_id}/{flowcell_id}/pod5'
        output:
            ubam = config.base_dir + config.ont_ubam + '/{sample_id}/{sample_id}_{flowcell_id}.bam',
            basecalling_qc_file = config.base_dir + config.qc_dir + '/basecalling_qc/{flowcell_id}_{sample_id}_cramino_qc.txt'
        resources:
            runtime=1440, mem_mb=150000, gpu=4, gpu_model='a100', disk_mb=50000
        threads: 30
        params:
            model = config.dorado_model,
            outdir = config.base_dir + config.ont_ubam + '/{sample_id}/'
        shell:
            """
            module load dorado/1.2.0
            module load nanopack/20231214
            scripts/basecalling.sh \
                --infile {input.pod5} \
                --outfile {output.ubam} \
                --outdir {params.outdir} \
                --model {params.model} \
                --qc-file {output.basecalling_qc_file} \
                --threads {threads}
            """
if not SKIP_BASECALL:
    TRIMMING_INPUT = rules.basecall.output.ubam
else:
    TRIMMING_INPUT = config.base_dir + config.ont_ubam + '/{sample_id}/{sample_id}_{flowcell_id}.bam'

rule trimming:
    input:  
        ubam = TRIMMING_INPUT
    output:  
        fastq = config.base_dir + config.pychopper_dir + '/{sample_id}/{sample_id}_{flowcell_id}.trimmed.fastq',
        trimming_qc_file = config.base_dir + config.qc_dir + '/trimming_qc/{sample_id}/{sample_id}_{flowcell_id}.tsv',
        read_stats = config.base_dir + config.qc_dir + '/trimming_qc/{sample_id}/{sample_id}_{flowcell_id}_read_stats.tsv',
        cramino_qc_file = config.base_dir + config.qc_dir + '/trimming_qc/cramino_stats/{sample_id}_{flowcell_id}_cramino_qc.txt'
    params:
        outdir = config.base_dir + config.pychopper_dir,
        kit = config.pychopper_kit,
        min_q_score = config.minimum_quality_score,
        min_length = config.minimum_read_length
    threads: 120
    resources:
        runtime=1440, mem_mb=120000, disk_mb=50000
    singularity:
        "./lrrna_1.2.sif"
    shell: 
        """
        scripts/trimming.sh \
            --infile {input.ubam} \
            --outfile {output.fastq}  \
            --kit {params.kit} \
            --min-quality {params.min_q_score} \
            --min-length {params.min_length} \
            --outdir {params.outdir} \
            --qc-file {output.trimming_qc_file} \
            --read-stats {output.read_stats} \
            --cramino-qc-file {output.cramino_qc_file} \
            --threads {threads} 
        """

if USE_FASTQ_SCREEN:
    rule fastq_screen:
        input:
            fastq = config.base_dir + config.pychopper_dir + '/{sample_id}/{sample_id}_{flowcell_id}.trimmed.fastq'
        output:
            txt = config.base_dir + config.qc_dir + '/fastq_screen/{sample_id}/{sample_id}_{flowcell_id}.trimmed_screen.txt',
            html = config.base_dir + config.qc_dir + '/fastq_screen/{sample_id}/{sample_id}_{flowcell_id}.trimmed_screen.html'
        params:
            outdir = config.base_dir + config.qc_dir + '/fastq_screen/{sample_id}',
            conf = config.fastq_screen_conf,
            subset = config.get('fastq_screen_subset', 100000)
        threads: 120
        resources:
            runtime=120, mem_mb=120000, disk_mb=50000
        singularity:
            "./lrrna_1.2.sif"
        shell:
            """
            scripts/fastq_screen.sh \
                --conf {params.conf} \
                --subset {params.subset} \
                --outdir {params.outdir} \
                --threads {threads} \
                --infile {input.fastq}
            """

rule alignment:
    input:  
        fastq = config.base_dir + config.pychopper_dir + '/{sample_id}/{sample_id}_{flowcell_id}.trimmed.fastq'
    output:  
        mapped_bam = config.base_dir + config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human_mapped.sorted.bam',
        #QC output files
        sirv_stats = config.base_dir + config.qc_dir + '/mapping_qc/{sample_id}/{sample_id}_{flowcell_id}_sirv_stats.txt' if USE_SIRV else [],
        mosdepth_sirv_summary = config.base_dir + config.qc_dir + '/mapping_qc/mosdepth_sirv/{sample_id}_{flowcell_id}.mosdepth.summary.txt' if USE_SIRV else [],
        cramino_sirv_stats = config.base_dir + config.qc_dir + '/mapping_qc/cramino_sirv_stats/{sample_id}_{flowcell_id}_cramino_qc.txt' if USE_SIRV else [],
        human_stats = config.base_dir + config.qc_dir + '/mapping_qc/{sample_id}/{sample_id}_{flowcell_id}_human_stats.txt',
        cramino_human_stats = config.base_dir + config.qc_dir + '/mapping_qc/cramino_human_stats/{sample_id}_{flowcell_id}_cramino_qc.txt',
        mosdepth_human_summary = config.base_dir + config.qc_dir + '/mapping_qc/mosdepth_human/{sample_id}_{flowcell_id}.mosdepth.summary.txt'
    threads: 120
    resources:
        runtime=900, mem_mb=200000, disk_mb=100000
    params:
        use_sirv = USE_SIRV,
        mapped_dir = config.base_dir + config.mapping_dir + '/{sample_id}/',
        mapped_qc_file = config.base_dir + config.qc_dir + '/mapping_qc/{sample_id}/{sample_id}_{flowcell_id}_stats.txt',
        human_fasta = config.genome,
        sirv_fasta = config.get('sirvome', ''),
        sirv_flag = lambda wildcards: f"--sirvome {config.get('sirvome', '')}" if USE_SIRV else "--skip-sirv"
    singularity:
        "./lrrna_1.2.sif"
    shell: 
        """
        scripts/alignment.sh \
            --infile {input.fastq} \
            --outfile {output.mapped_bam} \
            --outdir {params.mapped_dir} \
            --qc-file {params.mapped_qc_file} \
            --cramino-sirv-qc-file {output.cramino_sirv_stats} \
            --cramino-human-qc-file {output.cramino_human_stats} \
            --mosdepth-sirv-prefix {output.mosdepth_sirv_summary} \
            --mosdepth-human-prefix {output.mosdepth_human_summary} \
            --genome {params.human_fasta} \
            --threads {threads} \
            {params.sirv_flag} 
            
        """

rule stringtie:
    input:
        mapped_bam = config.base_dir + config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human_mapped.sorted.bam'
    output:
        sirv_stringtie_gtf = config.base_dir + config.stringtie_dir + '/{sample_id}/sirv/{sample_id}_{flowcell_id}_sirv.stringtie.gtf' if USE_SIRV else [],
        human_stringtie_gtf = config.base_dir + config.stringtie_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human.stringtie.gtf'
    threads: 30
    params:
        use_sirv = USE_SIRV,
        human_ref_gtf = config.human_genedb,
        sirv_ref_gtf = config.get('sirv_genedb', ''),
        sirv_stringtie_dir = config.base_dir + config.stringtie_dir + '/{sample_id}/sirv/',
        human_stringtie_dir = config.base_dir + config.stringtie_dir + '/{sample_id}/',
        assembly_mode = ASSEMBLY_MODE
    resources:
        runtime=120, mem_mb=12000, disk_mb=50000, slurm_partition='norm'
    singularity:
        "./lrrna_1.2.sif"
    shell: 
        """
        """ + ("scripts/sirv_stringtie_assembly.sh \
            --infile {input.mapped_bam} \
            --outfile {output.sirv_stringtie_gtf} \
            --outdir {params.sirv_stringtie_dir} \
            --ref_gtf {params.sirv_ref_gtf} \
            --threads {threads}\n" if USE_SIRV else "") + """
        scripts/human_stringtie_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.human_stringtie_dir} \
            --outfile {output.human_stringtie_gtf} \
            --ref_gtf {params.human_ref_gtf} \
            --mode {params.assembly_mode} \
            --threads {threads}
        """


rule isoquant:
    input:  
        mapped_bam = config.base_dir + config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human_mapped.sorted.bam'
    output:
        human_isoquant_gtf = config.base_dir + config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.transcript_models.gtf',
        human_isoquant_counts = config.base_dir + config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.discovered_transcript_counts.tsv',
    threads: 120
    params:
        use_sirv = USE_SIRV,
        prefix = '{sample_id}_{flowcell_id}',
        sirv_prefix = '{sample_id}_{flowcell_id}_sirv',
        human_fasta = config.genome,
        sirv_fasta = config.get('sirvome', ''),
        human_ref_gtf = config.human_genedb,
        sirv_ref_gtf = config.get('sirv_genedb', ''),
        sirv_isoquant_dir = config.base_dir + config.isoquant_dir + '/{sample_id}/sirv/',
        human_isoquant_dir = config.base_dir + config.isoquant_dir + '/{sample_id}/',
        assembly_mode = ASSEMBLY_MODE
    resources:
        runtime=4320, mem_mb=120000, disk_mb=50000
    singularity:
        "./lrrna_1.2.sif"
    shell: 
        """
        """ + ("scripts/sirv_isoquant_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.sirv_isoquant_dir} \
            --sirvome {params.sirv_fasta} \
            --genedb {params.sirv_ref_gtf} \
            --prefix {params.sirv_prefix} \
            --threads {threads}\n" if USE_SIRV else "") + """
        scripts/human_isoquant_assembly.sh \
            --infile {input.mapped_bam} \
            --outdir {params.human_isoquant_dir} \
            --genome {params.human_fasta} \
            --genedb {params.human_ref_gtf} \
            --prefix {params.prefix} \
            --mode {params.assembly_mode} \
            --threads {threads}
        """

rule merge:
    input:
        human_isoquant_gtf = config.base_dir + config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.transcript_models.gtf',
        human_isoquant_counts = config.base_dir + config.isoquant_dir + '/{sample_id}/{sample_id}_{flowcell_id}/{sample_id}_{flowcell_id}.discovered_transcript_counts.tsv',
        human_stringtie_gtf = config.base_dir + config.stringtie_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human.stringtie.gtf',
        mapped_bam = config.base_dir + config.mapping_dir + '/{sample_id}/{sample_id}_{flowcell_id}_human_mapped.sorted.bam'
    threads: 120
    params:
        prefix = '{sample_id}_{flowcell_id}',
        human_fasta = config.genome,
        human_ref_gtf = config.human_genedb,
        merge_dir = config.base_dir + config.merge_dir + '/{sample_id}',
        assembly_mode = ASSEMBLY_MODE,
        scripts_dir = config.script_dir,
    resources:
        runtime=480, mem_mb=400000, disk_mb=100000
    output:
        annotated_gtf = config.base_dir + config.merge_dir + '/{sample_id}/{sample_id}_{flowcell_id}.annotated.gtf'
    singularity:
        "./lrrna_1.2.sif"
    shell: 
        """
        scripts/transcript_merge.sh \
            --input_bam {input.mapped_bam} \
            --isoquant_gtf {input.human_isoquant_gtf} \
            --isoquant_counts {input.human_isoquant_counts} \
            --stringtie_gtf {input.human_stringtie_gtf} \
            --merge_dir {params.merge_dir} \
            --genome {params.human_fasta} \
            --ref_gtf {params.human_ref_gtf} \
            --prefix {params.prefix} \
            --outfile {output.annotated_gtf} \
            --assembly_mode {params.assembly_mode} \
            --scripts_dir {params.scripts_dir} \
            --threads {threads}
        """

rule create_multiqc_names:
    output:
        tsv = config.base_dir + config.qc_dir + '/multiqc_sample_names.tsv'
    run:
        import pandas as pd
        with open(output.tsv, 'w') as f:
            for idx, row in samples.iterrows():
                sample_id = row[0]
                flowcell_id = row[1]
                # creating a TSV file to rename sample names during report generation
                f.write(f"{sample_id}_{flowcell_id}_human_stats\t{sample_id}_human\n")
                if USE_SIRV:
                    f.write(f"{sample_id}_{flowcell_id}_sirv_stats\t{sample_id}_sirv\n")
                f.write(f"{sample_id}_{flowcell_id}\t{sample_id}\n")

rule multiqc_report:
    input:
        stats = expand(
            config.base_dir + config.qc_dir + '/mapping_qc/{sample_id}/{sample_id}_{flowcell_id}_human_stats.txt',
            zip,
            sample_id = samples.iloc[:,0].tolist(),
            flowcell_id = samples.iloc[:,1].tolist()
        ),
        fastq_screen = expand(
            config.base_dir + config.qc_dir + '/fastq_screen/{sample_id}/{sample_id}_{flowcell_id}.trimmed_screen.txt',
            zip,
            sample_id = samples.iloc[:,0].tolist(),
            flowcell_id = samples.iloc[:,1].tolist()
        ) if USE_FASTQ_SCREEN else [],
        name = config.base_dir + config.qc_dir + '/multiqc_sample_names.tsv'
    output:
        report = config.base_dir + config.qc_dir + '/multiqc_report.html'
    params:
        search_dir = config.base_dir + config.qc_dir,
        outdir = config.base_dir + config.qc_dir,
        multiqc_yaml = "config/multiqc_config.yaml"
    threads: 2
    resources:
        runtime=30, mem_mb=12000, disk_mb=5000
    singularity:
        "./lrrna_1.2.sif"
    shell:
        """
        multiqc {params.search_dir} \
            -o {params.outdir} \
            -n multiqc_report.html \
            --replace-names {input.name} \
            --config {params.multiqc_yaml} \
            --force \
        """
