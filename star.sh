#!/bin/bash

# Exit on error
set -euo pipefail

STAR_EXEC="/home/aghavidel/STAR/bin/Linux_x86_64/STAR"

# HARD CODED!!!
DATA_DIR="/home/aghavidel/data"
# Paths to the genome FA/GTF files
GENOME_FA="${DATA_DIR}/human_genome/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="${DATA_DIR}/human_genome/gencode.v45.annotation.gtf"
# STAR genome index directory
GENOME_INDEX="${DATA_DIR}/human_genome/index"
# Path to the FASTQ file we want to parse
READ_FASTQ="${DATA_DIR}/SRR16101430/SRR16101430_1.fastq"
# Read output directory
READ_OUTPUT_DIR="/home/aghavidel/read_outputs/"

THREADS="$1"

# Make the output dir if does not exist
mkdir -p "$READ_OUTPUT_DIR"

"${STAR_EXEC}" --runThreadN "$THREADS" \
     --genomeDir "${GENOME_INDEX}" \
     --readFilesIn "${READ_FASTQ}" \
     --outFileNamePrefix "${READ_OUTPUT_DIR}" \
     --outFilterMultimapNmax 20 \
     --outSAMtype SAM
