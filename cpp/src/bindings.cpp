#include "nucleotide_encoding.hpp"
#include "pam_search.hpp"
#include "tst.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_ternary_search_tree, m)
{

    m.doc() = "CRISPRitz C++ API bindings (pybind11)";

    // -------------------------------------------------------------------------
    // Ternary search tree API
    // -------------------------------------------------------------------------
    m.def("build_tree", &crispritz::build_tree, py::arg("sequence"), py::arg("pam"),
          py::arg("pam_length"), py::arg("pam_size"), py::arg("upstream"),
          py::arg("num_threads") = 1,
          "Build the TST index from an input sequence set using the provided PAM");
}
