#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pam_search.hpp"
#include "nucleotide_encoding.hpp"
#include "tst.hpp"


namespace py = pybind11;

PYBIND11_MODULE(_ternary_search_tree, m) {
    m.doc() = "CRISPRitz C++ API bindings (pybind11)";

    // -------------------------------------------------------------------------
    // Ternary search tree API 
    // -------------------------------------------------------------------------
    m.def("build_tree", &crispritz::build_tree, py::arg("name"), "Build a placeholder tree (stub).");


}
