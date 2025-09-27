#!/bin/bash

# Exit on error
set -euo pipefail

# ==== CONFIGURATION ====
THREADS=32

# Input files and directories
GENOME_FA="/mnt/data/aghavidel/human_genome/GRCh38.primary_assembly.genome.fa"
GTF="/mnt/data/aghavidel/human_genome/gencode.v45.annotation.gtf"


#ALZ
#FASTQ="/mnt/data/aghavidel/alz/SRR16101430_1.fastq"
FASTQ="/mnt/data/aghavidel/alz/SRR16101435_1.fastq"

#Healthy-Human
#FASTQ="/mnt/data/aghavidel/healthy_human/hg002_gm24385.mrna.R1.fastq"


GENOME_DIR="/mnt/data/aghavidel/human_genome/index"
OUT_PREFIX="/mnt/data/aghavidel/output/Alz"

# ==== Create index directory ====
mkdir -p "$GENOME_DIR"

# ==== Step 1: Generate genome index ====
#echo "Generating genome index..."
#/home/aghavidel/STAR/source/STAR  --runThreadN "$THREADS" \
#     --runMode genomeGenerate \
#     --genomeDir "$GENOME_DIR" \
#     --genomeFastaFiles "$GENOME_FA" \
#     --genomeSAindexNbases 14

echo "Current time: $(date)"

# ==== Step 2: Run STAR alignment ====
echo "Running STAR alignment... thread=8"
mkdir -p "$(dirname "$OUT_PREFIX")"
time /home/aghavidel/STAR/source/STAR --runThreadN "$THREADS" \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn $FASTQ \
     --outFileNamePrefix "$OUT_PREFIX" \
     --outFilterMultimapNmax 20 \
     --outSAMtype SAM

#echo "STAR run complete."

echo "Current time: $(date)"
