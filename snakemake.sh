#!/usr/bin/env bash
#SBATCH --time 10-0:00:00


module purge
module load apptainer
module load snakemake/7.32.4 

# Clone the biowulf snakemake profile
if [[ ! -d snakemake_profile ]]; then
    git clone https://github.com/NIH-HPC/snakemake_profile.git
fi

# Pull the containers
apptainer pull --disable-cache lrrna_latest.sif oras://quay.io/wellerca/lrrna

# Loading singularity
module load singularity/4.2.2

# Bind external directories on Biowulf
. /usr/local/current/singularity/app_conf/sing_binds

#updating permissions on bash scripts
chmod 777 scripts/trimming.sh
chmod 777 scripts/alignment.sh
chmod 777 scripts/human_isoquant_assembly.sh
chmod 777 scripts/sirv_isoquant_assembly.sh
chmod 777 scripts/human_stringtie_assembly.sh
chmod 777 scripts/sirv_stringtie_assembly.sh
chmod 777 scripts/transcript_merge.sh

# Run snakemake
# snakemake --profile snakemake_profile --use-singularity $@ -n
snakemake --profile snakemake_profile --use-singularity