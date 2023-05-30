#!/bin/bash

########################################################################
## Single-cell script to launch single-cell pipeline
##
## using: sbatch /mnt/beegfs/userdata/m_graa/pl_sc/script_projet_long/run_pipeline_analyse_int.sh
########################################################################
#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=longq
source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/m_graa/conda_env/single_cell2
module load singularity
#module load snakemake
path_to_configfile="/mnt/beegfs/userdata/m_graa/pl_sc/script_projet_long/parametres_analyses_integrees.yaml"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell/1.3"
snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile} #--unlock
conda deactivate