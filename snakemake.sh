#!/usr/bin/env bash
#SBATCH --time 02:00:00

################################################################################
# Snakemake Pipeline Launcher
# Explicity select the profile matching your system (Biowulf, Other HPCs or Local)
# Usage: ./snakemake.sh [biowulf|slurm|default]
################################################################################

# 1. Profile Selection
# User must explicitly provide the profile argument.
if [[ "$1" == "biowulf" || "$1" == "slurm" || "$1" == "default" ]]; then
    PROFILE="$1"
    shift  # Remove first argument so that "$@" contains only snakemake args
else
    echo "Error: You must specify an execution profile."
    echo "Usage: ./snakemake.sh <biowulf|slurm|default> [snakemake options]"
    echo ""
    echo "Available profiles:"
    echo "  biowulf : For NIH Biowulf cluster (loads modules automatically)"
    echo "  slurm   : For generic SLURM clusters"
    echo "  default : For local execution"
    exit 1
fi

################################################################################
# 2. Environment Loading
################################################################################

case "$PROFILE" in
    biowulf)
        # Load required modules for Biowulf
        module purge
        module load singularity/4.3.7 snakemake/7.32.4
        
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
if [[ ! -f "lrrna_1.2.sif" ]]; then
    echo "Container image not found. Downloading..."
    singularity pull oras://quay.io/datatecnica/lrrna:1.2
fi

################################################################################
# 3. Permissions & Execution
################################################################################

# Make scripts executable (using +x is safer than 777)
chmod +x scripts/*.sh

echo "Starting Snakemake pipeline with profile: $PROFILE"

# Run Snakemake
# We load default config

snakemake \
    --profile ./snakemake_profiles/$PROFILE \
    --config profile=$PROFILE \
    "$@"

# Exit with the same code as Snakemake (0 = success, 1 = error)
exit $?