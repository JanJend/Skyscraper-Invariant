cmake_minimum_required(VERSION 3.15)
project(skyscraper VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Need AIDA installed as a library
find_library(AIDA_LIBRARY NAMES aida
    PATHS ${AIDA_DIR}/build /usr/local/lib
    NO_DEFAULT_PATH
)
find_library(AIDA_LIBRARY NAMES aida)
if(NOT AIDA_LIBRARY)
    message(FATAL_ERROR "AIDA library not found. Build and install AIDA first.")
endif()
# Find Boost
find_package(Boost REQUIRED COMPONENTS timer chrono system)
#Find hdf5
find_package(HDF5 REQUIRED COMPONENTS CXX)
# Find CGAL
find_package(CGAL REQUIRED)

# Set build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set optimization flags for Release
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# Optimize file_reader.cpp even in Debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#    if(MSVC)
#        set_source_files_properties(src/file_reader.cpp PROPERTIES 
#            COMPILE_FLAGS "/O2")
#    else()
#        set_source_files_properties(src/file_reader.cpp PROPERTIES 
#            COMPILE_FLAGS "-O2")
#    endif()
# endif()
# Add _debug suffix to executables in Debug mode
set(CMAKE_DEBUG_POSTFIX "_debug")

# Change these paths to wherevery AIDA and Persistence-Algebra are.
set(AIDA_DIR "${CMAKE_SOURCE_DIR}/../AIDA" CACHE PATH "Path to AIDA directory")
set(PERSISTENCE_ALGEBRA_DIR "${CMAKE_SOURCE_DIR}/../Persistence-Algebra/include" CACHE PATH "Path to Persistence-Algebra directory")

# Verify that external dependencies exist
if(NOT EXISTS "${AIDA_DIR}/include/aida_interface.hpp")
    message(WARNING "AIDA not found at ${AIDA_DIR}. Set AIDA_DIR to correct path.")
endif()

if(NOT EXISTS "${PERSISTENCE_ALGEBRA_DIR}/include/grlina/graded_linalg.hpp")
    message(WARNING "Persistence-Algebra not found at ${PERSISTENCE_ALGEBRA_DIR}. Set PERSISTENCE_ALGEBRA_DIR to correct path.")
endif()

# Include directories
include_directories(include)
include_directories(${PERSISTENCE_ALGEBRA_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${AIDA_DIR}/include)
include_directories(${AIDA_DIR}/src)

# After project() declaration
option(BUILD_EXP "Build only hnf_at_origin and large_induced_indecomposables" OFF)

if(NOT BUILD_EXP)

    # Main HNF executable
    add_executable(hnf_main
        src/hnf.cpp
        src/uni_b1.cpp 
        src/hnf_at.cpp
        hnf_main.cpp
    )

    target_link_libraries(hnf_main ${AIDA_LIBRARY} ${Boost_TIMER_LIBRARY}
        ${Boost_CHRONO_LIBRARY}
        ${Boost_SYSTEM_LIBRARY} ${HDF5_LIBRARIES} CGAL::CGAL)
    set_target_properties(hnf_main PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

    # Makes filtered landscapes from .sky files
    add_executable(filt_landscape_from_sky
        src/filt_landscape.cpp
        src/file_reader.cpp
        filt_landscape_from_sky.cpp
    )
    set_target_properties(filt_landscape_from_sky PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
    target_link_libraries(filt_landscape_from_sky ${Boost_LIBRARIES})

    # Makes the module to quiver representation conversion
    add_executable(pres_to_quiver
        pres_to_quiver.cpp
    )
    set_target_properties(pres_to_quiver PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

    # Make the arrangement test
    add_executable(arrangement_test
        src/subdivision.cpp
        src/hnf_at.cpp
        src/uni_b1.cpp 
        arrangement_test.cpp
    )
    target_link_libraries(arrangement_test ${AIDA_LIBRARY} ${Boost_LIBRARIES} ${HDF5_LIBRARIES} CGAL::CGAL)
    set_target_properties(arrangement_test PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
endif()

add_executable(large_induced_indecomposables
    large_induced_indecomposables.cpp
)

target_link_libraries(large_induced_indecomposables ${AIDA_LIBRARY} ${Boost_TIMER_LIBRARY}
    ${Boost_CHRONO_LIBRARY}
    ${Boost_SYSTEM_LIBRARY})
set_target_properties(large_induced_indecomposables PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

add_executable(hnf_at_origin
    src/hnf_at.cpp
    src/uni_b1.cpp 
    hnf_at_origin.cpp
)
target_link_libraries(hnf_at_origin  ${Boost_TIMER_LIBRARY}
    ${Boost_CHRONO_LIBRARY}
    ${Boost_SYSTEM_LIBRARY} )
set_target_properties(hnf_at_origin PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

if(BUILD_EXP)
    install(TARGETS large_induced_indecomposables hnf_at_origin
        RUNTIME DESTINATION bin
    )
else()
    install(TARGETS hnf_main filt_landscape_from_sky pres_to_quiver arrangement_test large_induced_indecomposables hnf_at_origin
        RUNTIME DESTINATION bin
    )
endif()

