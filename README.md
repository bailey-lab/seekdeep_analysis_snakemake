# seekdeep_analysis_snakemake

This repository is attempting to be a general purpose seekdeep analysis
pipeline. Outputs include:
 - total read counts for every sample as a table
 - total read counts for every sample as a heatmap
 - complexity of infection (number of haplotypes per sample) as a table
 - complexity of infection (number of haplotypes per sample) as a heatmap
 - if applicable, relative abundance of amino acid point mutations as a table
 - if applicable, relative abundance of amino acid point mutations as a heatmap

## Installation:
### Install conda:
https://github.com/conda-forge/miniforge#mambaforge

### Create a conda environment and install snakemake there:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

### Setup your environment:
 - Change directory to a folder where you want to run the analysis
 - Download the analyze_seekdeep.smk file into this folder (e.g. with git clone)
 - Download the analyze_seekdeep.yaml file into same folder (e.g. git clone)

## Usage:
 - Edit the config.yaml file using the instructions in the comments. Use a text 
 editor that outputs unix line endings (e.g. vscode, notepad++, gedit, micro,
 emacs, vim, vi, etc.)
 - If snakemake is not your active conda environment, activate snakemake with:
```bash
conda activate snakemake
```
 - Run snakemake with:
```bash
snakemake -s analyze_seekdeep.smk --cores [your_desired_core_count]
```

