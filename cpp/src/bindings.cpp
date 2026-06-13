#include "nucleotide_encoding.hpp"
#include "pam_search.hpp"
#include "tst.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;

PYBIND11_MODULE(_ternary_search_tree, m) {
  m.doc() = "CRISPRitz C++ API bindings (pybind11)";

  // Map C++ runtime errors to Python RuntimeError so callers get a
  // descriptive exception instead of a hard crash.
  py::register_exception<std::runtime_error>(m, "TSTBuildError");
  py::register_exception<std::invalid_argument>(m, "TSTSearchError");

  // =========================================================================
  // TST index construction
  // =========================================================================
  m.def("build_tree", &crispritz::build_tree, py::arg("sequence"),
        py::arg("chr_name"), py::arg("pam_seq"), py::arg("pam_length"),
        py::arg("pam_limit"), py::arg("upstream"), py::arg("outdir"),
        py::arg("max_bulges") = 0, py::arg("num_threads") = 1,
        R"doc(
Build a Ternary Search Tree index for a single genomic sequence.

Parameters
----------
sequence : str
    Full genomic sequence (single chromosome, uppercase IUPAC).
chr_name : str
    Chromosome / contig identifier used in output filename(s).
pam_seq : str
    PAM-only string (e.g. ``"NGG"``), without guide placeholder Ns.
pam_length : int
    Total length of the PAM+guide pattern
    (e.g. 23 for ``NNNNNNNNNNNNNNNNNNNNNGG``).
pam_limit : int
    Length of the PAM portion only (e.g. 3 for ``NGG``).
upstream : bool
    True when the PAM precedes the guide (PAM-upstream, e.g. Cas12a ``TTTN``).
outdir : str
    Path to the directory where the genome index will be stored.
max_bulges : int, optional
    Maximum number of bulges; extra bases are extracted per site to support
    bulge-aware off-target search. Default 0.
num_threads : int, optional
    Number of OpenMP threads used during PAM search. Default 1.

Raises
------
TSTBuildError
    If no valid PAM sites are found or an output file cannot be written.
)doc");
}