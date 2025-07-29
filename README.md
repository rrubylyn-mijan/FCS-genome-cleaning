# FCS Genome Cleaning Pipeline

This repository contains detailed, step-by-step commands used to screen, clean, mask, and assess contaminant adapter sequences in wheat genome assemblies using:

- FCS-adaptor and FCS-GX tools from NCBI
- Singularity containers
- BEDTools, SAMtools, BLAST+, and EMBOSS

## Instructions

All commands are listed in the file:  
**`FCS_genome_cleaning_pipeline.sh`**

The steps include:
- Screening with `run_fcsadaptor.sh`
- Cleaning with `fcs.py`
- Masking contaminants
- Comparing overlap with gene annotations
- Translating affected sequences
- Running local BLAST and analyzing alignments

---

## Maintainer

This guide was written by **Ruby Mijan**

- Tools developed by: [NCBI FCS Team](https://github.com/ncbi/fcs)
- This repository only organizes and documents how to use FCS tools effectively in an HPC environment using Singularity.

> This is not the official FCS repository. This is a practical guide to simplify the execution of NCBI’s FCS-adaptor and FCS-GX tools to screen and clean wheat genome assemblies.