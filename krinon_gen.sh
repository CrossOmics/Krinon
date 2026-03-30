#!/bin/bash

# Exit on error
set -euo pipefail

ALIGNER_ROOT="$(dirname "$(readlink -f "${BASH_SOURCE}")")"
ALIGNER_EXEC="${ALIGNER_ROOT}/build/RNAAlignRefactored"

# HARD CODED!!!
DATA_DIR="/home/aghavidel/data"
# Paths to the genome FA/GTF files
GENOME_FA="${DATA_DIR}/human_genome/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="${DATA_DIR}/human_genome/gencode.v45.annotation.gtf"
# Krinon genome index directory
GENOME_INDEX="${DATA_DIR}/human_genome/refactored_index"

# Make the output dir if does not exist
mkdir -p "$GENOME_INDEX"

# Genome generattion
"${ALIGNER_EXEC}" \
     --mode GenomeGenerate \
     --buildMethod SAIS \
     --genomeFile "${GENOME_FA}" \
     --gtfFile "${GENOME_GTF}" \
     --genomeGenerateFileStoreDir "${GENOME_INDEX}"
