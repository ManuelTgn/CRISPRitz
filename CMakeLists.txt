set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
cmake_minimum_required(VERSION 3.19)
project(crispritz LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find pybind11 (provided by pip/conda)
find_package(pybind11 CONFIG REQUIRED)

# Build python extension module: crispritz/_crispritz.* (import as crispritz._crispritz)
pybind11_add_module(_ternary_search_tree
    cpp/src/bindings.cpp
    cpp/src/tst.cpp
    cpp/src/pam_search.cpp
    cpp/src/nucleotide_encoding.cpp
)

target_include_directories(_ternary_search_tree PRIVATE
    cpp/include
)

# Optional: OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_ternary_search_tree PRIVATE OpenMP::OpenMP_CXX)
endif()

# Absolute path to the Python package directory
set(PY_PACKAGE_DIR "${CMAKE_SOURCE_DIR}/src/crispritz")

install(TARGETS _ternary_search_tree
    LIBRARY DESTINATION ${PY_PACKAGE_DIR} 
)
