#!/bin/bash

# Exit on error
set -euo pipefail

# HARD CODED!!!
DATA_DIR="/ssd1/arvin/Aligner"
STAR_EXEC="/home/arvin/Aligner/STAR/source/STAR"
# STAR_EXEC="/home/arvin/Aligner/STAR/bin/Linux_x86_64/STAR"

GENOME_DIR="${DATA_DIR}/genomes"
OUT_DIR="${DATA_DIR}/outputs"
READ_DIR="${DATA_DIR}/reads"

GENOME_FA="${GENOME_DIR}/human_genome/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="${GENOME_DIR}/human_genome/gencode.v45.annotation.gtf"
GENOME_INDEX="${GENOME_DIR}/human_genome/index"

# FASTQ="/mnt/data/aghavidel/alz/SRR16101435_1.fastq"
READ_FASTQ="${READ_DIR}/SRR16101430_1.fastq"

THREADS="$1"

echo "Running STAR alignment... thread=${THREADS}"
time "${STAR_EXEC}" --runThreadN "$THREADS" \
     --genomeDir "$GENOME_INDEX" \
     --readFilesIn $READ_FASTQ \
     --outFileNamePrefix "$OUT_DIR" \
     --outFilterMultimapNmax 20 \
     --outSAMtype SAM

echo "STAR run complete."

echo "Current time: $(date)"
