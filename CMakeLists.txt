cmake_minimum_required(VERSION 3.19)
project(crispritz LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Export compile commands for IDEs and static analysis tools
# ------------------------------------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # -std=c++17, not -std=gnu++17

# ------------------------------------------------------------------------------
# Default build type
# Scikit-build-core sets this from the pip/cmake CLI; default to Release so a
# plain cmake .. & make also produces an optimised build
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type specified - defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
find_package(pybind11 CONFIG REQUIRED)

# OpenMP is required: the PAM search and TST construction are parallelised with 
# #pragma omp and will silently fall back to single-threaded execution on 
# MacOS/clang without it, but we want a clear error rather than a silent 
# correctness surprise
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------------------------
# Python extension module
# ------------------------------------------------------------------------------
pybind11_add_module(_ternary_search_tree
    cpp/src/bindings.cpp
    cpp/src/tst.cpp
    cpp/src/tst_utils.cpp
    cpp/src/pam_search.cpp
    cpp/src/nucleotide_encoding.cpp
)

target_include_directories(_ternary_search_tree PRIVATE
    cpp/include
)

target_link_libraries(_ternary_search_tree PRIVATE
    OpenMP::OpenMP_CXX
)

# ------------------------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------------------------
target_compile_options(_ternary_search_tree PRIVATE
    # Common flags for all compilers
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter   # many pybind11 shims trigger this
    >
 
    # Release-only: full optimisation + auto-vectorisation hints
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang,AppleClang>,$<CONFIG:Release>>:
        -O3
        -march=native           # SIMD (AVX2/SSE4) where the build host supports it
        -funroll-loops
        -ffast-math             # safe for genomic integer-heavy code; no NaN/Inf
    >
 
    # Debug-only: sanitisers catch the kind of UB common in index arithmetic
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang,AppleClang>,$<CONFIG:Debug>>:
        -O0
        -g3
        -fsanitize=address,undefined
    >
 
    # MSVC equivalents (Windows CI)
    $<$<CXX_COMPILER_ID:MSVC>:
        /W4
        /O2
    >
)
 
# The sanitiser flags must also appear at link time
target_link_options(_ternary_search_tree PRIVATE
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang,AppleClang>,$<CONFIG:Debug>>:
        -fsanitize=address,undefined
    >
)

# ------------------------------------------------------------------------------
# Install
# ------------------------------------------------------------------------------
set(PY_PACKAGE_DIR "${CMAKE_SOURCE_DIR}/src/crispritz")

install(TARGETS _ternary_search_tree
    LIBRARY DESTINATION ${PY_PACKAGE_DIR} 
)
