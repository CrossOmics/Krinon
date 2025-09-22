### RNA Aligner

```
--mode
GenomeGenerate, ReadAlign, both

--genomeFile
path to genome fasta file

--readFile
path to read file (fastq)

--readFile2
path to read file 2 (fastq) [for paired-end]

--readType
paired, single

--gtfFile
test.gtf

--kMerSize
length of k-mer in GenomeIndex, default 14

--outPutDir
output directory

--sjdbOverhang
same as STAR

--threads
1

--outFilterMultimapMax
the maximum number of multiple alignments allowed for each read, default 10

--genomeGenerateFileStoreDir
directory to store genome index
```
