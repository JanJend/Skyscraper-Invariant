cmake_minimum_required(VERSION 3.15)
project(skyscraper VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Find Boost
find_package(Boost REQUIRED COMPONENTS timer system)

# Set build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set optimization flags for Release
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# Add _debug suffix to executables in Debug mode
set(CMAKE_DEBUG_POSTFIX "_debug")

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

# Main HNF executable
add_executable(hnf_main
    hnf_main.cpp
    src/filt_landscape.cpp
)

target_link_libraries(hnf_main ${Boost_LIBRARIES})
set_target_properties(hnf_main PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

# Makes filtered landscapes from .sky files
add_executable(filt_landscape_from_sky
    filt_landscape_from_sky.cpp
    src/filt_landscape.cpp
)
set_target_properties(filt_landscape_from_sky PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

# Makes the module to quiver representation conversion
add_executable(pres_to_quiver
    pres_to_quiver.cpp
)
set_target_properties(pres_to_quiver PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

target_link_libraries(filt_landscape_from_sky ${Boost_LIBRARIES})