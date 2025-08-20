cmake_minimum_required(VERSION 3.22)
project(RNAAlignRefactored)

set(CMAKE_CXX_STANDARD 17)

add_executable(RNAAlignRefactored   main.cpp
                                    genome/Genome.cpp
                                    genome/GenomeIndex.cpp
                                    genome/GenomeIndexPrefix.cpp
                                    genome/SuffixArray.cpp
                                    genome/SuffixAutomatonKMerBuilder.cpp
                                    genome/SJDB.cpp
                                    readAlign/ReadAligner.cpp
                                    readAlign/SeedMapping.cpp
                                    readAlign/StitchingManagement.cpp
                                    genome/WindowScanningKMerBuilder.cpp
                                    utils/PackedArray.cpp
                                    )
