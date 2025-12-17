cmake_minimum_required(VERSION 3.15)
project(CRISPRitz VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find OpenMP
find_package(OpenMP)

# Library source files
set(PAM_LIB_SOURCES
    cpp/src/nucleotide_encoding.cpp
    cpp/src/pam_search.cpp
)

# Library header files
set(PAM_LIB_HEADERS
    cpp/include/crispritz/nucleotide_encoding.hpp
    cpp/include/crispritz/pam_search.hpp
)

# Create library
add_library(pam_lib ${PAM_LIB_SOURCES} ${PAM_LIB_HEADERS})

# Set include directories
target_include_directories(pam_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(pam_lib PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found - parallel search enabled")
else()
    message(WARNING "OpenMP not found - search will run single-threaded")
endif()

# Compiler-specific optimization flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building with optimizations enabled")
    if(MSVC)
        target_compile_options(pam_lib PRIVATE /O2 /W4)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.14)
            target_compile_options(pam_lib PRIVATE /arch:AVX2)
        endif()
    else()
        target_compile_options(pam_lib PRIVATE 
            -O3 
            -march=native 
            -mtune=native
            -Wall 
            -Wextra
            -Wpedantic
        )
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building with debug symbols")
    if(MSVC)
        target_compile_options(pam_lib PRIVATE /W4 /Zi)
    else()
        target_compile_options(pam_lib PRIVATE -g -O0 -Wall -Wextra)
    endif()
endif()

# ============================================================================
# Testing
# ============================================================================

enable_testing()

# Test executable
add_executable(test_pam_search tests/test_basic.cpp)
target_link_libraries(test_pam_search PRIVATE pam_lib)

# Add test to CTest
add_test(NAME BasicPAMSearchTests COMMAND test_pam_search)

# Set test properties
set_tests_properties(BasicPAMSearchTests PROPERTIES
    TIMEOUT 60
    LABELS "unit"
)

# ============================================================================
# Optional: Example executable (commented out since no main yet)
# ============================================================================

# Uncomment when main.cpp is ready:
# if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/examples/main.cpp)
#     add_executable(pam_generator examples/main.cpp)
#     target_link_libraries(pam_generator PRIVATE pam_lib)
# endif()

# ============================================================================
# Installation rules
# ============================================================================

install(TARGETS pam_lib
    EXPORT PAMSearchTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(FILES ${PAM_LIB_HEADERS}
    DESTINATION include
)

# Install CMake config files for find_package support
install(EXPORT PAMSearchTargets
    FILE PAMSearchTargets.cmake
    NAMESPACE PAMSearch::
    DESTINATION lib/cmake/PAMSearch
)

# ============================================================================
# Print configuration summary
# ============================================================================

message(STATUS "")
message(STATUS "=== PAM Search Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "OpenMP support: ${OpenMP_CXX_FOUND}")
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "===============================")
message(STATUS "")