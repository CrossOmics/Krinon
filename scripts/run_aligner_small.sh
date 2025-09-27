#!/bin/bash

# HARD CODED!!!
DATA_DIR="/ssd1/arvin/Aligner"

ALIGNER_ROOT="$(dirname "$(readlink -f "${BASH_SOURCE}")")/.."

GENOME_DIR="${DATA_DIR}/genomes"
OUT_DIR="${DATA_DIR}/outputs"
READ_DIR="${DATA_DIR}/reads"

ALIGNER_EXEC="${ALIGNER_ROOT}/build/RNAAlignRefactored"

GENOME_FA="${GENOME_DIR}/human_genome/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="${GENOME_DIR}/human_genome/gencode.v45.annotation.gtf"
GENOME_INDEX="${GENOME_DIR}/human_genome/new_index"

READ_FASTQ="${READ_DIR}/small.fastq"

THREADS="$1"

echo "Current time: $(date)"

echo "Running alignment... thread=${THREADS}"
mkdir -p "$(dirname "$OUT_PREFIX")"
time "${ALIGNER_EXEC}" --threads "${THREADS}" \
     --mode ReadAlign \
     --genomeGenerateFileStoreDir "${GENOME_INDEX}" \
     --readFile "${READ_FASTQ}" \
     --outPutDir "${OUT_DIR}" \
     --outFilterMultimapMax 20 

echo "code run complete."

echo "Current time: $(date)"
