# [Arvin Ghavidel - 3/12/2026]

## Current State
- I have verified that the genome index building appears fine, and built the index for human genome.
- I have also confirmed that I am no longer seeing the significant memory spike per thread. It appears that it actually might have been a bug due to invalid casting, casuing some of the loaded sizes to be much bigger than they should have been.
- There is a problem with the logic for making the output SAM file larger works. Currently, the file starts small with 1e9 bytes, and
is then expanded when needed. Originally, I explicitly tried to avoid doing something like this since I wanted to avoid using any lock for writing to
the output SAM file or changing the file size, but I can see why you would prefer to change that and until we benchmark it, I will avoid suggesting you to change this.
With that said, the logic is faulty. The `BUSERR` happens because of the incorrect way `mremap` is called. There are two aspects to a memory mapped file:
    - The virtual adderss space in the calling process which we write data to
    - The acutal file on disk, which the OS guarantees to mirror the address space unless the process aborts or we call `munmap`

    One cannot resize one of these without resizing the other. The `BUSERR` happens because currently, `ensureSize` only resizes the virtual memory address space and does not update the file on disk. As such, the moment we pass the original boundry of the address space, the OS reports that we have crossed into an unmapped region of the disk (note that this is _NOT_ a `SEGFAULT`, since the memory address itself is valid), and thus aborts with a `BUSERR`.
    
    The correct way to fix this is to resize _BOTH_ the address space with `mremap` and the file on disk with `fallocate`. The trouble here is that `fallocate` is quite slow, but since the number of calls to it will be logarithmic in the input size, it hopefully amortizes over the whole runtime.

## Changes
- Some changes to how file names are handled (e.g. see `main.cpp`). Instead of manually joining paths, I am using (and strictly recommend) using `std::filesystem`, since we are using C++17 anyway.
- Changes to `ensureSize` in `MemoryMappedFile.hpp` to fix the `BUSERR`.

## Suggestions
- When loading the genome, I see that there are a lot of cases where we want to return some sort of error code. Consider documenting these and
using an `enum` instead. Check the `TODO` in `GenomeIndex_build.cpp`.

# [Arvin Ghavidel - 3/15/2026]
- Added link time optimization to the CMake build file just to squeeze a bit more out of compiler optimization. Helps a little bit (single-threaded output
went from 40M/hour to around 53M/hour).

# [Arvin Ghavidel - 3/21/2026]
- Merged changes from `main` by Yuxuan. I have verified that the single-threaded performance is back to its expected value (around 200 million per hour on
the SPHERE server), but the multi-threaded performance is very bad. As such, I have concluded that the current writing method to the output is not efficient.
Currently, the IO works as follows:
    - Aligner threads are spawned, and grab hold of a shared pointer on the `mmap`ed input file.
    - Threads race with each other on grabbing the next read, parse it, increment the pointer and then align. To prevent race on the pointer, the 
    `readFileLock_` mutex is used to guard `loadFromFastq`, putting any other thread to sleep while we are reading the mapped value.
    - After alignment is done by calling the `StitchingManagement` object, the threads create the string to output, then race to write to the output
    using `outputFileLock_` to prevent mangling.
This kills the performance, we are reaping very little benefit from having `mmap`ed the gigantic file, since we are not parallelizing over the input.
- On a server with NVMe, a single thread can saturate write throughput, thus we change the IO format to the following:
    - The file is mapped, and then chunked into `N` partitions. Each thread is positioned to start reading from the initial part of the partition (care
    will be taken to prevent mangling the read lines).
    - Each thread aligns lock-free (the chunks do not overlap), and writes its output to buffer.
    - When the output buffer of a thread is full, the result will be passed to a lock-free queue.
    - On the other end, a single thread consumes the buffer and writes the result, cleaning the buffer and returning it to the caller.

>**The OOM Fix**  
>For huge read files, we map the file in chunks (the maximum size of the chunk will be user defined). Each chunk will be processed normally, then the next chunk is mapped and we start again.

# [Arvin Ghavidel - 3/22/2026]
- Most of the changes above have been implemented. With some major changes to the previous proposed pipeline.

## New Aligner Pipeline

- The input SAM file is partitioned into chunks of pre-specified size (currently, 32 GB).
- Each partition is mapped into the aligner process exactly once, aligned and released. We then proceed partition-by-partition until we are done.
- Each aligner, loads a smaller region of the mapped partition into its input buffer, aligns it, and then pushes the result into its output buffer.
- When the output buffer is full, the result is submitted to a single therad dedicated to _only_ writing to the output file. The output file is no
longer mapped and lives on disk.

Some details:
- The current implementation uses a lock-free queue to communicate with the writer thread. In order to prevent contention and to make sure that the
write throughput is saturated with a single thread, the output buffer size needs to be quite big. A few MB (current 8 MB) seems to do the trick for
both NVMe and HDD storage.
- The input data is no longer chunked per-thread. Threads grab buffer-fulls of data from the current offset into the mapped file and then increment (this
is basically how Yuxuan originally implemented this). The problem with the chunked version (which I originally proposed under the `refine` branch last 
year) is that it suffers from stragglers, of which there seems to be always at least one, meaning that my initial intuition about the alignment procedure
being more uniform seems to be very incorrect. For this reason, I restored Yuxuan's original implementation for the read scanners.
- One difference with Yuxuan's implementaiton however, is that instead of using fixed length buffers, the buffers are now basically just strings which
have been preallocated before we start. This gives no performance penalty, while making sure that a `SEGFAULT` cannot happen in the rare even that we
get a buffer overflow for a read that has been matched with a huge amount of things. It also means that we don't have to manage these buffers directly.

## Changes
- A multi-producer, single-consumer, lock-free queue has been implemented in `io/MPSC.hpp`. This is quite bare-bones, but I didn't want to add an extra
dependency just for that. With benefit of hindsight, given the large buffers, even a potatoey lock would do just fine. The current implementation also
has a (not _fully_ lock-free) buffer pool for managing the output buffers.
- Some functions that previuosly _returned_ strings, now receive a reference to a string (which is actually a buffer). These include:
    - `convertToSAM` (`Transcript.h`)
    - `process` and `convertToResult` (`Stitching.h`)
- I had to rewrite `convertToSAM` in `Transcript.cpp` slightly to make it zero-copy (basically, just added an inline function for some string operation 
and replaced a bunch of string operations with `append`)
- All `const char*` buffers of any kind (especially in `Stitching`) have been removed. We don't need to bother doing these things manually.
- `ReadAligner` and `OutputSAM` have be rewritten to accomodate the pipeline above. They are much simpler things now.
- Replace some plain pointers with smart pointers out of sheer paranoia ...

## TODOs
- Some extra logic needs to be added for handling partitioned reads.
- After reading the code, I have concluded that we are not at all prepared to handle paired-end reads. The current (and very likely previous) implementations can easily
output garbage if used as is. Thus I think will be pretty difficult to solve, so I am postponing that.

# [Arvin Ghavidel - 3/25/2026]

## Changes
- The output SAM file creation now has `O_TRUNC` added. If it exists already, it will be truncated.
- The `MemoryMappedFile` class has been renamed to `MemeoryMappedInput` and all write-capable functionalities for it have been removed (we don't need them).
- I have ditched the partitioning strategy. We'll just map the entire file and let the kernel handle paging by itself. I have added `madvise` to the new
implementation to hint for sequential access.
