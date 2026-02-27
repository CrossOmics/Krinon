### Compile: 
```
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```
### Run:
```
--mode
GenomeGenerate,ReadAlign

--buildMethod
SAIS(Suggested),Traditional

--genomeFile
YOUR_GENOME_REF (xxx.fa)

--gtfFile
GTF_FILE (xxx.gtf)

--readFile
test_short1.fastq

--readType
single,paired(under development)

--kMerSize
14 (same as STAR)

--genomeGenerateFileStoreDir
./test/

--outPutDir
./test/

--readBufferSize
(ReadBufferSize of each threads,default: 5e7)
when there are many threads, it is recommended to set a bigger value

--outputBufferSize
(OutputBufferSize of each threads,default: 5e7)
same as readBufferSize

--threads
# of threads to use, default: 1
```
For more parameters, please check Parameter.cpp.

I have not conducted large-scale data testing on this version, so there might be some correctness or speed issues that I haven't identified. I appreciate your understanding. I will address them promptly once I am informed. :)