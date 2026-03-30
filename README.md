# Krinon is an RNA aligner producing almost the same output as STAR.

The package is still under development. So far, we only support one-pass alignment for bulk-RNA reads.

Some works are in our todo list now.

1. Two-pass
2. Chimeric
3. Single-cell read alignment.
4. HiFi long-read

## Quick Build and Run

To build Krinon:
```sh
mkdir build
cmake -S src/ -B build/ -DCMAKE_BUILD_TYPE=Release
pushd build
make -j
popd
```

Use the scripts `krinon_align` and `krinon_gen` to quickly test things. Paths to any files required
can be configured there.

```sh
# Generate genome index with FA/GTF files
./krinon_gen.sh
# Invoke the aligner with a FASTQ file and 32 thread
./krinon_align.sh 32
```
