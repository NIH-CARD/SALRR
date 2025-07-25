#!/bin/bash
# properties = {"type": "single", "rule": "basecall", "local": false, "input": ["rawdata/NABEC_KEN-1069_FTX_RNA/NABEC_KEN-1069_FTX_RNA/20250114_2319_2A_PAW72725_5119b730/pod5"], "output": ["rawdata/ONT_UBAM/NABEC_KEN-1069_FTX_RNA/NABEC_KEN-1069_FTX_RNA_20250114_2319_2A_PAW72725_5119b730.bam"], "wildcards": {"sample_id": "NABEC_KEN-1069_FTX_RNA", "flowcell_id": "20250114_2319_2A_PAW72725_5119b730"}, "params": {"model": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"}, "log": [], "threads": 1, "resources": {"mem_mb": 300000, "mem_mib": 286103, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "runtime": 4320, "gpu": 4, "gpu_model": "v100x"}, "jobid": 5, "cluster": {}}

# if lscratch exists use it for tempdir
if [[ -d "/lscratch/$SLURM_JOB_ID" ]] ; then
    tmp="/lscratch/$SLURM_JOB_ID/tmp"
    mkdir "$tmp"
    export TMPDIR="$tmp"
fi

cd /vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA && /usr/local/apps/snakemake/conda/envs/7.32.4/bin/python3.11 -m snakemake --snakefile '/vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA/snakefile' --target-jobs 'basecall:sample_id=NABEC_KEN-1069_FTX_RNA,flowcell_id=20250114_2319_2A_PAW72725_5119b730' --allowed-rules 'basecall' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=300000' 'mem_mib=286103' 'disk_mb=1000' 'disk_mib=954' 'gpu=4' --wait-for-files '/vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA/.snakemake/tmp.y_odsa8r' 'rawdata/NABEC_KEN-1069_FTX_RNA/NABEC_KEN-1069_FTX_RNA/20250114_2319_2A_PAW72725_5119b730/pod5' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'params' 'code' 'mtime' 'software-env' 'input' --skip-script-cleanup  --conda-frontend 'mamba' --singularity-args ' --cleanenv' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 240 --scheduler 'greedy' --scheduler-solver-path '/usr/local/apps/snakemake/conda/envs/7.32.4/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1

