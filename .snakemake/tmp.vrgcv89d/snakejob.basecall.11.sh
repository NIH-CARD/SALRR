#!/bin/bash
# properties = {"type": "single", "rule": "basecall", "local": false, "input": ["/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-00-34_FTX_RNA/NABEC_SH-00-34_FTX_RNA/20250114_2321_2E_PAY64520_2d2a8c92/pod5"], "output": ["/data/CARD_AUX/LRS_temp/NABEC_RNA/SNAKEMAKE_TEST/ONT_UBAM/NABEC_SH-00-34_FTX_RNA/NABEC_SH-00-34_FTX_RNA_20250114_2321_2E_PAY64520_2d2a8c92.bam"], "wildcards": {"sample_id": "NABEC_SH-00-34_FTX_RNA", "flowcell_id": "20250114_2321_2E_PAY64520_2d2a8c92"}, "params": {"model": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "outdir": "/data/CARD_AUX/LRS_temp/NABEC_RNA/SNAKEMAKE_TEST/ONT_UBAM/NABEC_SH-00-34_FTX_RNA/"}, "log": [], "threads": 1, "resources": {"mem_mb": 300000, "mem_mib": 286103, "disk_mb": 50000, "disk_mib": 47684, "tmpdir": "<TBD>", "runtime": 7200, "gpu": 4, "gpu_model": "v100x"}, "jobid": 11, "cluster": {}}

# if lscratch exists use it for tempdir
if [[ -d "/lscratch/$SLURM_JOB_ID" ]] ; then
    tmp="/lscratch/$SLURM_JOB_ID/tmp"
    mkdir "$tmp"
    export TMPDIR="$tmp"
fi

cd /vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA && /usr/local/apps/snakemake/conda/envs/7.32.4/bin/python3.11 -m snakemake --snakefile '/vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA/snakefile' --target-jobs 'basecall:sample_id=NABEC_SH-00-34_FTX_RNA,flowcell_id=20250114_2321_2E_PAY64520_2d2a8c92' --allowed-rules 'basecall' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=300000' 'mem_mib=286103' 'disk_mb=50000' 'disk_mib=47684' 'gpu=4' --wait-for-files '/vf/users/CARD_AUX/LRS_temp/NABEC_RNA/CARDlongread_ONT_long_read_RNA/.snakemake/tmp.vrgcv89d' '/data/CARD_AUX/LRS_temp/NABEC_RNA/NABEC_SH-00-34_FTX_RNA/NABEC_SH-00-34_FTX_RNA/20250114_2321_2E_PAY64520_2d2a8c92/pod5' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'params' 'input' 'mtime' 'code' --skip-script-cleanup  --conda-frontend 'mamba' --singularity-args ' --cleanenv' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 240 --scheduler 'greedy' --scheduler-solver-path '/usr/local/apps/snakemake/conda/envs/7.32.4/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1

