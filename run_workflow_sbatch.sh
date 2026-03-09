#!/usr/bin/bash
#SBATCH --job-name=spatial_multiome         # Job name
#SBATCH --mem=200G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=18:00:00                       # Time limit 12 hours
#SBATCH --output=./results/spatial_multiome_%j.log               # Standard output and error log
#SBATCH --cpus-per-task=16 					# num cores

source activate spatial_multiome
# Deeptools should really be included in the environment
module load deeptools

snakemake -s workflow.snakefile --cores 16 --rerun-incomplete
