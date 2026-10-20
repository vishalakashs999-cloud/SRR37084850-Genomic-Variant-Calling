#!/usr/bin/env bash
# SRR37084850 — E. coli genomic variant calling
# Curated from the original shell history; not independently rerun or validated.
# Requirements: paired trimmed FASTQ files, FastQC, BWA, SAMtools, GATK 4.6.2.0,
# Java, SnpEff with a configured custom Ecoliref database.
set -euo pipefail

WORKDIR="/mnt/d/genomic"
GATK="/home/vishal/gatk-4.6.2.0/gatk"
SNPEFF_DIR="/home/vishal/snpEff"
REF="GCF_000005845.2_ASM584v2_genomic.fna"
cd "$WORKDIR"

# Reference download and indexing
wget -nc "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/${REF}.gz"
gzip -dk "${REF}.gz"  # keep compressed reference copy
cp "$REF" referencegenome.fasta
bwa index referencegenome.fasta

# Quality control of pre-existing trimmed paired-end reads
fastqc SRR37084850_1_val_1.fq SRR37084850_2_val_2.fq

# Align reads
bwa mem referencegenome.fasta SRR37084850_1_val_1.fq SRR37084850_2_val_2.fq > demoaligned.sam
samtools flagstat demoaligned.sam > alignment_flagstat.txt

# SAM/BAM processing: name sort, fixmate, coordinate sort, duplicate removal
samtools view -b -o demoaligned.bam demoaligned.sam
samtools sort -n -o demoaligned_sort.bam demoaligned.bam
samtools fixmate -m demoaligned_sort.bam demoaligned_fixmate.bam
samtools sort -o demoaligned_sortfixmate.bam demoaligned_fixmate.bam
samtools markdup -r demoaligned_sortfixmate.bam demoalignedfinal.bam

# Add read groups and index BAM
python3 "$GATK" AddOrReplaceReadGroups \
  -I demoalignedfinal.bam -O demoaligned_final1.bam \
  -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 1
samtools index demoaligned_final1.bam
python3 "$GATK" CreateSequenceDictionary -R referencegenome.fasta
samtools faidx referencegenome.fasta
samtools coverage demoaligned_final1.bam > alignment_coverage.txt

# Original HaplotypeCaller invocation used -ERC GVCF.
# This produces a gVCF, not a conventional variant-only VCF.
python3 "$GATK" HaplotypeCaller \
  -R referencegenome.fasta -I demoaligned_final1.bam \
  -O demooutput.vcf -ERC GVCF

# IMPORTANT: The historical command filtered DP > 30 as "Fail".
# This direction is retained faithfully; it is NOT an endorsed threshold.
# Verify intended filtering and GVCF handling before reusing scientifically.
python3 "$GATK" VariantFiltration \
  -R referencegenome.fasta -V demooutput.vcf \
  -O FilteredVariants_new.vcf \
  --filter-expression 'DP > 30' --filter-name 'Fail'
python3 "$GATK" SelectVariants \
  -R referencegenome.fasta -V FilteredVariants_new.vcf \
  --exclude-filtered -O FilteredVariants_new3.vcf

# Annotation assumes custom SnpEff Ecoliref database is already built.
# Original history used 'java -jar snpEff.jar build -gtf22 ... -V Ecoliref'
# after preparing data/Ecoliref/sequences.fa and genes.gtf and editing snpEff.config.
cp FilteredVariants_new3.vcf "$SNPEFF_DIR/FilteredVariants_new3.vcf"
cd "$SNPEFF_DIR"
java -jar snpEff.jar Ecoliref FilteredVariants_new3.vcf > demoannotation.vcf
