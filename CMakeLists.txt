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
# plain cmake .. & make also produces an optimised build.
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
# macOS/clang without it, but we want a clear error rather than a silent
# correctness surprise.
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------------------------
# Compiler / linker flag lists
# These must be defined before any target that references them.
# ------------------------------------------------------------------------------
set(CRISPRITZ_COMPILE_FLAGS
    # Warnings – applied unconditionally for GCC / Clang / AppleClang
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

# Sanitiser flags must also be passed at link time
set(CRISPRITZ_LINK_FLAGS
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang,AppleClang>,$<CONFIG:Debug>>:
        -fsanitize=address,undefined
    >
)

# ------------------------------------------------------------------------------
# Shared object library for all C++ sources
# Compiling the sources once into a static library avoids redundant compilation
# across the extension module and every test executable.
# ------------------------------------------------------------------------------
add_library(crispritz_core STATIC
    cpp/src/nucleotide_encoding.cpp
    cpp/src/offtarget.cpp
    cpp/src/output_writer.cpp
    cpp/src/pam_search.cpp
    cpp/src/profile_data.cpp
    cpp/src/search_configuration.cpp
    cpp/src/tst.cpp
    cpp/src/tst_search.cpp
    cpp/src/tst_utils.cpp
)

target_include_directories(crispritz_core PUBLIC cpp/include)
target_link_libraries(crispritz_core PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(crispritz_core PRIVATE ${CRISPRITZ_COMPILE_FLAGS})
target_link_options(crispritz_core PRIVATE ${CRISPRITZ_LINK_FLAGS})

# <filesystem> requires an explicit link on GCC < 9 and some Clang versions.
# Linking stdc++fs is harmless on compilers that bundle it into libstdc++.
target_link_libraries(crispritz_core PUBLIC
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
)

# ------------------------------------------------------------------------------
# Python extension module
# ------------------------------------------------------------------------------
pybind11_add_module(_ternary_search_tree
    cpp/src/bindings.cpp
)

target_link_libraries(_ternary_search_tree PRIVATE crispritz_core)
target_compile_options(_ternary_search_tree PRIVATE ${CRISPRITZ_COMPILE_FLAGS})
target_link_options(_ternary_search_tree PRIVATE ${CRISPRITZ_LINK_FLAGS})

# ------------------------------------------------------------------------------
# C++ unit tests (only built when CTest is enabled)
# ------------------------------------------------------------------------------
include(CTest)  # defines BUILD_TESTING (ON by default)

if(BUILD_TESTING)

    # Helper macro – avoids repeating the same boilerplate for every test target.
    # Usage: add_crispritz_test(<target> <source_file>)
    macro(add_crispritz_test TARGET SOURCE)
        add_executable(${TARGET} ${SOURCE})
        target_link_libraries(${TARGET} PRIVATE crispritz_core)
        target_compile_options(${TARGET} PRIVATE ${CRISPRITZ_COMPILE_FLAGS})
        target_link_options(${TARGET} PRIVATE ${CRISPRITZ_LINK_FLAGS})
        # Register with CTest; run from the build directory so relative paths work
        add_test(
            NAME    ${TARGET}
            COMMAND ${TARGET}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
        # Treat a non-zero exit code as a failure (default CTest behaviour).
        set_tests_properties(${TARGET} PROPERTIES
            FAIL_REGULAR_EXPRESSION "\\[FAIL\\]"
        )
    endmacro()

    add_crispritz_test(test_nucleotide_encoding
        cpp/tests/test_nucleotide_encoding.cpp
    )

    add_crispritz_test(test_offtarget
        cpp/tests/test_offtarget.cpp
    )

    add_crispritz_test(test_output_writer
        cpp/tests/test_output_writer.cpp
    )

    add_crispritz_test(test_pam_search
        cpp/tests/test_pam_search.cpp
    )

    add_crispritz_test(test_profile_data
        cpp/tests/test_profile_data.cpp
    )

    add_crispritz_test(test_search_configuration
        cpp/tests/test_search_configuration.cpp
    )

    add_crispritz_test(test_tst_utils
        cpp/tests/test_tst_utils.cpp
    )

    add_crispritz_test(test_tst
        cpp/tests/test_tst.cpp
    )

    add_crispritz_test(test_tst_search
        cpp/tests/test_tst_search.cpp
    )

endif()  # BUILD_TESTING

# ------------------------------------------------------------------------------
# Install
# ------------------------------------------------------------------------------
set(PY_PACKAGE_DIR "${CMAKE_SOURCE_DIR}/src/crispritz")

install(TARGETS _ternary_search_tree
    LIBRARY DESTINATION ${PY_PACKAGE_DIR}
)
