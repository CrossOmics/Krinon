# cmake_minimum_required(VERSION 3.22)
cmake_minimum_required(VERSION 3.16)
project(RNAAlignRefactored)

# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)


add_executable(RNAAlignRefactored main.cpp
        #genome/Genome.cpp
        #genome/GenomeIndex.cpp
        #genome/SuffixArray.cpp
        #genome/SJDB.cpp
        #readAlign/ReadAligner.cpp
        #readAlign/SeedMapping.cpp
        #readAlign/Stitching.cpp
        #readAlign/ReadAlignMultiThread.cpp
        #readAlign/ReadFile.cpp
        #utils/PackedArray.cpp
        #utils/argparse.hpp
        #utils/Parameters.cpp
        #utils/PackedIndexArray.cpp
        #utils/MemoryMappedFile.hpp
        utilsRefactored/seqFunctions.cpp
        log/ErrorRecord.cpp
        log/TimeReport.cpp
        utilsRefactored/PackedArray.cpp
        utilsRefactored/timeFunctions.cpp
        io/Parameters.cpp
        io/MemoryMappedFile.hpp
        io/ReadScanner.cpp
        io/OutputSAM.cpp
        genomeRefactored/Genome.cpp
        genomeRefactored/GenomeIndex_build.cpp
        genomeRefactored/GenomeIndex_find.cpp
        genomeRefactored/SuffixArray.cpp
        genomeRefactored/SuffixArray_build_traditional.cpp
        genomeRefactored/SuffixArray_build_SAIS.cpp
        genomeRefactored/SuffixArray_build_CAPS.cpp
        genomeRefactored/SuffixArrayKMerMap.cpp
        genomeRefactored/SJDB.cpp
        readAlignRefactored/Stitching.cpp
        readAlignRefactored/SeedMapping.cpp
        readAlignRefactored/Transcript.cpp
        readAlignRefactored/ReadMerger.cpp
        readAlignRefactored/ReadAligner.cpp
)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RNAAlignRefactored PRIVATE -march=native)
elseif (MSVC)
    target_compile_options(RNAAlignRefactored PRIVATE /arch:AVX2)
endif ()

target_link_libraries(RNAAlignRefactored PUBLIC OpenMP::OpenMP_CXX)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(RNAAlignRefactored PRIVATE -Wall -Wextra -Wpedantic -Wno-unused-variable -Wno-unused-parameter -Werror)
endif ()