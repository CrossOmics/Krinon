cmake_minimum_required(VERSION 3.16)
project(RNAAlignRefactored)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -pthread)

add_subdirectory(3rd-party/plog)
add_subdirectory(3rd-party/yaml-cpp)
add_executable(RNAAlignRefactored   main.cpp
                                    genome/Genome.cpp
                                    genome/GenomeIndex.cpp
                                    genome/SuffixArray.cpp
                                    genome/SJDB.cpp
                                    readAlign/ReadAligner.cpp
                                    readAlign/SeedMapping.cpp
                                    readAlign/StitchingManagement.cpp
                                    readAlign/ReadAlignMultiThread.cpp
                                    readAlign/ReadFile.cpp
                                    utils/PackedArray.cpp
                                    utils/argparse.hpp
                                    utils/Parameters.cpp
                                    )
target_link_libraries(RNAAlignRefactored plog::plog yaml-cpp::yaml-cpp)
