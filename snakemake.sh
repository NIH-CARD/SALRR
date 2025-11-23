#!/usr/bin/env bash

################################################################################
# Snakemake Pipeline Launcher
# automatically detects environment (Biowulf vs. Generic HPC vs. Local)
# Usage: ./snakemake.sh [biowulf|slurm|default]
################################################################################

# 1. Auto-Detect Logic
# If user provided an argument (e.g., "biowulf"), use it. Otherwise use "auto".
PROFILE="${1:-auto}"

if [[ "$PROFILE" == "auto" ]]; then
    # Check specifically for Biowulf
    if [[ -n "$SLURM_CLUSTER_NAME" ]] && [[ "$SLURM_CLUSTER_NAME" == "biowulf" ]]; then
        PROFILE="biowulf"
        echo "Detected Biowulf environment - using 'biowulf' profile"
    
    # Check for generic SLURM (checking if 'sinfo' command exists)
    elif command -v sinfo &> /dev/null; then
        PROFILE="slurm"
        echo "Detected SLURM cluster - using 'slurm' profile"
    
    # Fallback to local execution
    else
        PROFILE="default"
        echo "No cluster detected - using 'default' local profile"
    fi
fi

################################################################################
# 2. Environment Loading
################################################################################

case "$PROFILE" in
    biowulf)
        # Load required modules for Biowulf
        module purge
        module load apptainer singularity/4.2.2 snakemake/7.32.4
        
        # Load Biowulf-specific singularity bindings if they exist
        if [ -f /usr/local/current/singularity/app_conf/sing_binds ]; then
             . /usr/local/current/singularity/app_conf/sing_binds
        fi
        ;;
        
    slurm)
        # For generic HPC, assume Conda, but try loading module if missing
        if ! command -v snakemake &> /dev/null; then
            echo "Please activate the environment: conda activate lrrna"
            exit 1
        fi 
        if ! command -v singularity &> /dev/null && ! command -v apptainer &> /dev/null; then
        echo "Error: singularity/apptainer not found. Please load module: module load singularity"
        exit 1
        fi
        ;;
        
    default)
        # For local execution, strictly require user to have environment active
        if ! command -v snakemake &> /dev/null; then
            echo "Error: snakemake not found. Please run: conda activate lrrna"
            exit 1
        fi
        ;;
esac

# Ensure container exists (Standard across all profiles)
if [[ ! -f "lrrna_0.9.sif" ]]; then
    echo "Container image not found. Downloading..."
    singularity pull oras://quay.io/datatecnica/lrrna:0.9
fi

################################################################################
# 3. Permissions & Execution
################################################################################

# Make scripts executable (using +x is safer than 777)
chmod +x scripts/*.sh

echo "Starting Snakemake pipeline with profile: $PROFILE"

# Run Snakemake
# We load default config AND resources config
snakemake \
    --profile ./snakemake_profiles/$PROFILE \
    --configfile config/default.yaml \
    --config environment=$PROFILE \
    "$@"

# Exit with the same code as Snakemake (0 = success, 1 = error)
exit $?