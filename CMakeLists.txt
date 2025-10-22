cmake_minimum_required(VERSION 3.15)
project(skyscraper VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Find Boost
find_package(Boost REQUIRED COMPONENTS timer system)

# Change these paths to wherevery AIDA and Persistence-Algebra are.
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
include_directories(${Boost_INCLUDE_DIRS})

# Main HNF executabl
add_executable(hnf_main
    hnf_main.cpp
    src/hnf.cpp
    src/filt_landscape.cpp
)

target_link_libraries(hnf_main ${Boost_LIBRARIES})

# Makes filtered landscapes from .sky files
add_executable(filt_landscape_from_sky
    filt_landscape_from_sky.cpp
    src/filt_landscape.cpp
)
target_link_libraries(filt_landscape_from_sky ${Boost_LIBRARIES})