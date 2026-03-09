#!/usr/bin/bash
#SBATCH --job-name=make_bowtie_ref         # Job name
#SBATCH --partition=pe2                      # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=64G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=6:00:00                       # Time limit 12 hours
#SBATCH --output=./make_bowtie_ref_%j.log               # Standard output and error log

source activate spatial_multiome

bowtie2-build Mus_musculus.GRCm39.dna.primary_assembly.fa.gz GrCm39
