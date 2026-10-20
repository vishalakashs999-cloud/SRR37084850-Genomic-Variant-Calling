# Genomic Variant Calling and Functional Annotation — SRR37084850

## Project Overview

A bioinformatics workflow for processing publicly available *Escherichia coli* sequencing data, identifying genomic variants, and annotating the resulting variant calls.

**SRA Accession:** SRR37084850  
**Reference Genome:** GCF_000005845.2 (E. coli K-12 MG1655)

## Analysis Workflow

1. **Quality Control — FastQC:** Assessed sequencing-read quality.
2. **Read Alignment — BWA-MEM:** Aligned paired-end sequencing reads to the reference genome.
3. **BAM Processing — SAMtools:** Converted SAM to BAM, performed sorting, mate-information processing, duplicate removal, and indexing.
4. **Variant Calling — GATK HaplotypeCaller:** Generated genomic variant calls in VCF format.
5. **Variant Filtering — GATK:** Applied variant filtering and selected variants for downstream analysis.
6. **Variant Annotation — SnpEff:** Used a custom reference database to annotate variants.

## Project Outputs

| File | Description |
|---|---|
| `FilteredVariants_new3.zip` | Archived filtered VCF output |
| `demoannotation.zip` | Archived SnpEff-annotated VCF output |
| `snpEff.config` | SnpEff configuration file |
| `snpEff_genes.txt` | Supporting gene-related output |

## Tools and Technologies

Linux (Ubuntu/WSL), FastQC, BWA-MEM, SAMtools, GATK, SnpEff, NCBI SRA, and NCBI RefSeq.

## Project Scope

This repository documents a practical genomic variant-calling workflow, including reference preparation, sequence alignment, BAM processing, variant calling, filtering, and annotation.

The archived VCF files contain the generated outputs. Variant counts, filtering performance, and biological interpretations should be verified directly from those files before reporting specific findings.
