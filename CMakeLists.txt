cmake_minimum_required(VERSION 3.15)
project(skyscraper VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# External dependency paths - can be overridden by user
set(AIDA_DIR "${CMAKE_SOURCE_DIR}/../AIDA" CACHE PATH "Path to AIDA directory")
set(PERSISTENCE_ALGEBRA_DIR "${CMAKE_SOURCE_DIR}/../Persistence-Algebra" CACHE PATH "Path to Persistence-Algebra directory")

# Verify that external dependencies exist
if(NOT EXISTS "${AIDA_DIR}/aida_interface.hpp")
    message(WARNING "AIDA not found at ${AIDA_DIR}. Set AIDA_DIR to correct path.")
endif()

if(NOT EXISTS "${PERSISTENCE_ALGEBRA_DIR}/include/grlina/graded_linalg.hpp")
    message(WARNING "Persistence-Algebra not found at ${PERSISTENCE_ALGEBRA_DIR}. Set PERSISTENCE_ALGEBRA_DIR to correct path.")
endif()

# Include directories
include_directories(include)
include_directories(${AIDA_DIR})
include_directories(${PERSISTENCE_ALGEBRA_DIR}/include)

# Main HNF executable
add_executable(hnf_main
    hnf_main.cpp
    src/hnf.cpp
    src/filt_landscape.cpp
)

# Landscape from sky executable (doesn't need hnf)
add_executable(filt_landscape_from_sky
    filt_landscape_from_sky.cpp
    src/filt_landscape.cpp
)