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
