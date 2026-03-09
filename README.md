
# SPACE-Tag-preprocessing
Snakemake workflow for preprocessing of SPACE-Tag raw data for downstream analysis

## Introduction
This pipeline processes FASTQ files from a SPACE-Tag experiment. In this experiment, R1 contains the spatial barcode (16 bp) and Unique Molecular Identifier (UMI, 12 bp), while R2 contains the genomic sequences. The workflow employs various tools for trimming, mapping, and deduplicating the sequencing reads. Subsequently, the mapped reads are used for peak calling and feature counting. The final output is a set of files for downstream analysis such as using Seurat.

## Table of Contents

1. [Installation](#installation)
2. [Setting up a run](#setting-up-a-run)
3. [Configuring workflow options](#configuring-workflow-options)
4. [Running the workflow](#running-the-workflow)

## Installation

First, clone the main branch of the repository to your local machine.

```
git clone https://github.com/nygctech/SPACE-Tag-preprocessing.
cd SPACE-Tag-preprocessing
```

The requirements to run the workflow can be found in `environment.yml`. Create a conda environment with the required dependencies and activate it.

```
conda env create -f environment.yml
conda activate SPACE-Tag
```

You're ready to start!

## Setting up a run

Input files and options to the snakemake pipeline are specified in the `config.yaml` config file.

#### `sample_fastqs` - Input data

```
sample_fastqs:
  sample_name:
    # Read 1 (16bp Visium barcode + 12bm UMI)
    R1: path/to/read1.fastq.gz
    # Read 2 (gDNA only)
    R2: path/to/read2.fastq.gz
```

#### `sample_alignment_jsons` - Tissue image alignment

```
sample_alignment_jsons:
  # Path to SpaceRanger fiducial alignment json - if you want to do automated alignment (currently not implemented), put 'none' 
  sample_name: path/to/alignment.json
```

The `json` file specifies which spots are part of the tissue and is used to create the filtered feature matrix files in the directory `filtered_peak_bc_matrix`. The pipeline expects an input json file in the 10X SpaceRanger format (see [10X Genomics documentation](https://www.10xgenomics.com/support/software/space-ranger/analysis/inputs/image-fiducial-alignment)). 

#### `run_type` - Run type specification

This can be set to either `processing` or `debug`. If set to `processing`, `SPACE-Tag` will run the preprocessing workflow and generate output objects for downstream analyiss; if set to `debug`, `SPACE-Tag` will instead run a debug workflow and generate plots for QC and debugging. Currently, only pA/pT artifact filtering (see below) is implemented in the `debug` workflow. 

```
### Whether to do full run for processing or just debug/QC
# two options: 'processing', 'debug'
run_type: 'debug'
```

#### `ref` - Reference genome

`SPACE-Tag` uses the `bowtie2` aligner and expects references in a compatible format. If you have a `genome.fa` or `genome.fa.gz` file that contains the genome you would like to align to, run the following command to generate a `bowtie2`-compatible reference:

```
bowtie2-build genome.fa {prefix}
```

Replace `prefix` with a suitable name. See the script `refs/GrCm39/make_bowtie_ref.sh` for an example of this. References are specified in the config file using:

```
# Path and prefix for bowtie2 genome reference
ref: path/to/ref/and/prefix
```

## Configuring workflow options

The workflow has a number of options that can be adjusted. These options are specified in the `config.yaml` configuration file, where each option is specified as a `param:value` pair. Options are explained below:

#### `mode` - Paired-end or Single-end mode

Whether to use paired-end or single-end for alignment. If using the SPACE-Tag protocol, this should be set to `se` by default. If set to `pe`, a compatible `r1_format` must also be set (see below). 

```
### single-end or paired-end mode
# options: 'se', 'pe'
mode: 'se'
```

#### `r1_format` - Structure of Read 1

Structure of read one, containing spatial barcode, UMI (optional) and genomic read (optional). This is specified in the `umi_tools` regex format (see [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method)). If using the SPACE-Tag protocol, this should be set to `'(?P<cell_1>.{16})(?P<umi_1>.{12})'` by default, which specifies a 28bp read 1 containing a 16bp spatial barcode followed by a 12bp UMI.

```
# structure of read 1
# if using SPACE-Tag, default is '(?P<cell_1>.{16})(?P<umi_1>.{12})'
r1_format: '(?P<cell_1>.{16})(?P<umi_1>.{12})'
```


#### `umi-dedup`

```
# Whether to use UMIs for deduplication
umi_dedup: False
```

Boolean value that specifies whether deduplication should be performed by UMI or by genomic alignment position.

### Peak calling options

#### `macs2_fdr`

```
### Options for MACS2 peak calling
# FDR threshold for peak calling
macs2_fdr: 0.01
```

False discovery rate for MACS2 peak calling. Reasonable options are `0.01`, `0.05`, etc.

#### `macs2_genomesize`

```
# Effective genome size
# Defaults for species:
#   hs: 2.7e9
#   mm: 1.87e9
#   ce: 9e7
#   dm: 1.2e8
macs2_genomesize: 1.87e9
```

Effective genome size for MACS2 peak calling. Default sizes for common species are listed above.

### Additional configuration options

#### `remove_artifacts`

The current spatial CUT&Tag chemistry is vulnerable to spurious peaks at genomic regions containing A/T homopolymers. If this option is set to `True`, peaks with genomic poly A/T sequences above a certain length within a certain distance of the peak will be filtered out. Thresholds can be set below. If you are uncertain which parameters to use for filtering, set `run_type = 'debug'` above; this will generate a plot showing the number of peaks filtered out at a range of threshold values, as well as a plot showing the degree of nucleotide bias which can help in determining if artifactual peaks are present. 

```
# Remove artifact pA/pT peaks
remove_artifacts: True
```

#### `whitelist`

```
# Path for text file containing barcode whitelist
whitelist: refs/visium-v1.txt
```

The whitelist text file contains all valid spatial barcodes, with one barcode on each line. For the Visium v1 protocol, this file is provided in this repository.

#### `spot_coords`

```
# Path for text file containing spatial row/column indices of valid barcodes
spot_coords: refs/visium-v1_coordinates.txt
```

The spot coordinates file is a tab-delimited file containing:

Column 1: all valid spatial barcodes (these must match the barcodes specified in `whitelist` above) \
Column 2: 1-based column index of the spatial barcode  \
Column 3: 1-based row index of the spatial barcode 

### Artifactual peak filtering options

If `remove_artifacts` is set to `True`, these options set the thresholds below which artifactual peaks will be filtered out. If there is a poly A/T sequence of length greater than `poly_a_t_threshold` within `window_length` of a peak, it will be filtered out.

#### `window_length`

The size in bp of the window upstream and downstream of peaks to look for poly A/T genomic sequences. 

```
window_length: 70
```

#### `poly_a_t_threshold`

The length of poly A/T genomic sequence above which a peak will be filtered out. 

```
poly_a_t_threshold: 10
```

#### `genome_fasta`

If `remove_artifacts` is set to `True` to enable filtering out artifactual peaks at genomic poly A/T regions, you must also specify the path to a `genome.fa` or `genome.fa.gz` file. This should be the same genome version as was used to create the `bowtie2` reference to ensure correct peak filtering.

```
# This needs to be specified if remove_artifacts is set!
genome_fasta: path/to/genome/fasta
```

## Running the workflow

Once the `config.yaml` file has been properly configured, run the following command to start the workflow, replacing `{num_cores}` with an integer number of threads to use for execution:

```
snakemake -s workflow.snakefile -j {num_cores}
```

Alternatively, if you are running this in a SLURM-managed HPC system, you can use the `run_workflow_sbatch.sh` script to run the pipeline in the background using sbatch. In that case, run:

```
sbatch run_workflow_sbatch.sh
```

This should work assuming you have set up the environment correctly; if not, check the sbatch parameters at the top of the script.

(Optional) before running the full workflow, you can run snakemake with options `-np` to execute a dry run:

```
snakemake -s workflow.snakefile -np
```
