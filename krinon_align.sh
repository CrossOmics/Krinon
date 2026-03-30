#!/bin/bash

# Exit on error
set -euo pipefail

ALIGNER_ROOT="$(dirname "$(readlink -f "${BASH_SOURCE}")")"
ALIGNER_EXEC="${ALIGNER_ROOT}/build/RNAAlignRefactored"

# HARD CODED!!!
DATA_DIR="/home/aghavidel/data"
# Krinon genome index directory
GENOME_INDEX="${DATA_DIR}/human_genome/refactored_index"
# Path to the FASTQ file we want to parse
READ_FASTQ="${DATA_DIR}/SRR16101430/SRR16101430_1.fastq"
# Read output directory
READ_OUTPUT_DIR="/home/aghavidel/read_outputs/"

THREADS="$1"

# Make the output dir if does not exist
mkdir -p "$READ_OUTPUT_DIR"

# Alignment
"${ALIGNER_EXEC}" --threads "${THREADS}" \
     --mode ReadAlign \
     --readFile "${READ_FASTQ}" \
     --genomeGenerateFileStoreDir "${GENOME_INDEX}" \
     --outPutDir "${READ_OUTPUT_DIR}"
