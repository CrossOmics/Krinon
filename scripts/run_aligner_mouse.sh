
THREADS=32


GENOME_FA="/mnt/data/aghavidel/mouse_genome/GRCm38_genomic.fna"
GTF="/mnt/data/aghavidel/mouse_genome/GRCm38_genomic.gtf"


GENOME_DIR="/mnt/data/aghavidel/mouse_genome/new_index"
OUT_PREFIX="/mnt/data/aghavidel/new_output/mouse_erv/"


echo "Current time $(date)"

echo "Generating genome index..."
/home/aghavidel/RNA-Aligner-main/src/build/RNAAlignRefactored  --threads "$THREADS" \
     --mode GenomeGenerate \
     --genomeGenerateFileStoreDir "$GENOME_DIR" \
     --genomeFile "$GENOME_FA" \
      --kMerSize 13


echo "Genome generation complete"

echo "Current time: $(date)"

# ==== Step 2: Run alignment ====
#echo "Running alignment... thread=8"
#mkdir -p "$(dirname "$OUT_PREFIX")"
#time /home/aghavidel/RNA-Aligner-main/src/build/RNAAlignRefactored --threads "$THREADS" \
#     --mode ReadAlign\
#     --genomeGenerateFileStoreDir "$GENOME_DIR" \
#     --readFile $FASTQ \
#     --outPutDir "$OUT_PREFIX" \
#     --outFilterMultimapMax 20 

#echo "code run complete."

echo "Current time: $(date)"
