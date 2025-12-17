cmake_minimum_required(VERSION 3.22)
project(RNAAlignRefactored)

set(CMAKE_CXX_STANDARD 17)

add_executable(RNAAlignRefactored   main.cpp
                                    genome/Genome.cpp
                                    genome/GenomeIndex.cpp
                                    genome/SuffixArray.cpp
                                    genome/SJDB.cpp
                                    readAlign/ReadAligner.cpp
                                    readAlign/SeedMapping.cpp
                                    readAlign/Stitching.cpp
                                    readAlign/ReadAlignMultiThread.cpp
                                    readAlign/ReadFile.cpp
                                    utils/PackedArray.cpp
                                    utils/argparse.hpp
                                    utils/Parameters.cpp
                                    utils/PackedIndexArray.cpp
                                    utils/MemoryMappedFile.hpp
                                    utilsRefactored/seqFunctions.cpp
                                    log/ErrorRecord.cpp
                                    utilsRefactored/timeFunctions.hpp
                                    utilsRefactored/PackedArray.cpp
                                    )
