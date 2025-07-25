#!/usr/bin/env bash
#SBATCH --time 2-0:00:00


module purge
module load snakemake/7.32.4 

# Clone the biowulf snakemake profile
if [[ ! -d snakemake_profile ]]; then
    git clone https://github.com/NIH-HPC/snakemake_profile.git
fi

# Run snakemake
# snakemake --profile snakemake_profile $@ -n
snakemake --profile snakemake_profile -j 2