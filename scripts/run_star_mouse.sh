#!/bin/bash

# Exit on error
set -euo pipefail

# ==== CONFIGURATION ====
THREADS=32

# Input files and directories
GENOME_FA="/mnt/data/aghavidel/mouse_genome/GRCm38_genomic.fna"
GTF="/mnt/data/aghavidel/mouse_genome/GRCm38_genomic.gtf"

GENOME_DIR="/mnt/data/aghavidel/mouse_genome/index"
OUT_PREFIX="/mnt/data/aghavidel/output/mouse_ibd/"

echo "Current time: $(date)"

# ==== Create index directory ====
mkdir -p "$GENOME_DIR"

# ==== Step 1: Generate genome index ====
echo "Generating genome index..."
/home/aghavidel/STAR/source/STAR  --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FA" \
     --genomeSAindexNbases 14


echo "Current time: $(date)"

# ==== Step 2: Run STAR alignment ====
#echo "Running STAR alignment... thread=8"
#mkdir -p "$(dirname "$OUT_PREFIX")"
#time /home/aghavidel/STAR/source/STAR --runThreadN "$THREADS" \
#     --genomeDir "$GENOME_DIR" \
#     --readFilesIn $FASTQ \
#     --outFileNamePrefix "$OUT_PREFIX" \
#     --outFilterMultimapNmax 20 \
#     --outSAMtype SAM

#echo "STAR run complete."

echo "Current time: $(date)"
