#include "nucleotide_encoding.hpp"
#include "pam_search.hpp"
#include "search.hpp"
#include "tst.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;

PYBIND11_MODULE(_ternary_search_tree, m)
{
    m.doc() = "CRISPRitz C++ API bindings (pybind11)";

    // Map C++ runtime errors to Python RuntimeError so callers get a
    // descriptive exception instead of a hard crash.
    py::register_exception<std::runtime_error>(m, "TSTBuildError");
    py::register_exception<std::invalid_argument>(m, "TSTSearchError");

    // =========================================================================
    // TST index construction
    // =========================================================================
    m.def("build_tree", &crispritz::build_tree, py::arg("sequence"), py::arg("chr_name"),
          py::arg("pam_seq"), py::arg("pam_length"), py::arg("pam_limit"), py::arg("upstream"),
          py::arg("outdir"), py::arg("max_bulges") = 0, py::arg("num_threads") = 1,
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

    // =========================================================================
    // OffTarget result type
    // =========================================================================
    py::class_<crispritz::OffTarget>(m, "OffTarget",
                                     R"doc(
A single CRISPR off-target hit.

Attributes
----------
bulge_type : str
    One of ``"X"`` (no bulge), ``"DNA"``, ``"RNA"``, or ``"DNA,RNA"``.
guide_seq : str
    Guide RNA sequence with PAM ``N``-placeholders appended (or prepended).
target_seq : str
    Genomic target sequence; mismatch bases are lower-case.
chromosome : str
    Chromosome / contig name.
genomic_pos : int
    1-based genomic start position.
cluster_pos : int
    Cluster-normalised position (accounts for bulge length offsets).
strand : str
    ``'+'`` (forward) or ``'-'`` (reverse).
mismatches : int
    Number of substitution mismatches.
bulge_size : int
    Total number of bulge bases (DNA + RNA).
total_score : int
    ``mismatches + bulge_size``.
)doc")
        .def_readonly("bulge_type", &crispritz::OffTarget::bulge_type)
        .def_readonly("guide_seq", &crispritz::OffTarget::guide_seq)
        .def_readonly("target_seq", &crispritz::OffTarget::target_seq)
        .def_readonly("chromosome", &crispritz::OffTarget::chromosome)
        .def_readonly("genomic_pos", &crispritz::OffTarget::genomic_pos)
        .def_readonly("cluster_pos", &crispritz::OffTarget::cluster_pos)
        .def_readonly("strand", &crispritz::OffTarget::strand)
        .def_readonly("mismatches", &crispritz::OffTarget::mismatches)
        .def_readonly("bulge_size", &crispritz::OffTarget::bulge_size)
        .def_readonly("total_score", &crispritz::OffTarget::total_score)
        .def("__repr__",
             [](const crispritz::OffTarget& ot)
             {
                 return "<OffTarget chr=" + ot.chromosome +
                        " pos=" + std::to_string(ot.genomic_pos) + " strand=" + ot.strand +
                        " mm=" + std::to_string(ot.mismatches) +
                        " bulge=" + std::to_string(ot.bulge_size) + " type=" + ot.bulge_type + ">";
             });

    // =========================================================================
    // PositionalProfile result type
    // =========================================================================
    py::class_<crispritz::PositionalProfile>(m, "PositionalProfile",
                                             R"doc(
Per-guide positional mismatch and bulge counts.

All positional vectors have length ``guide_length``.
Histogram vectors have length ``max_mismatches + 1``.

Attributes
----------
guide : str
    Guide sequence with PAM N-placeholder appended (or prepended).
mm_per_pos : list[int]
    Total mismatch count at each guide position.
dna_per_pos : list[int]
    Total DNA-bulge count at each guide position.
rna_per_pos : list[int]
    Total RNA-bulge count at each guide position.
mm_counts : list[int]
    Number of off-target hits with exactly k mismatches (0 bulges),
    for k = 0 .. max_mismatches.
dna_counts : list[list[int]]
    dna_counts[k][j] = hits with k mismatches and j+1 DNA bulges.
rna_counts : list[list[int]]
    rna_counts[k][j] = hits with k mismatches and j+1 RNA bulges.
joint_counts : list[list[list[int]]]
    joint_counts[k][bd][br] = hits with k mismatches, bd DNA bulges,
    and br RNA bulges.
nuc_per_pos : list[list[list[int]]]
    nuc_per_pos[k][n][i] = count of nucleotide n (0=A,1=C,2=G,3=T) at
    position i in targets with k mismatches.
)doc")
        .def_readonly("guide", &crispritz::PositionalProfile::guide)
        .def_readonly("mm_per_pos", &crispritz::PositionalProfile::mm_per_pos)
        .def_readonly("dna_per_pos", &crispritz::PositionalProfile::dna_per_pos)
        .def_readonly("rna_per_pos", &crispritz::PositionalProfile::rna_per_pos)
        .def_readonly("mm_counts", &crispritz::PositionalProfile::mm_counts)
        .def_readonly("dna_counts", &crispritz::PositionalProfile::dna_counts)
        .def_readonly("rna_counts", &crispritz::PositionalProfile::rna_counts)
        .def_readonly("joint_counts", &crispritz::PositionalProfile::joint_counts)
        .def_readonly("nuc_per_pos", &crispritz::PositionalProfile::nuc_per_pos);

    // =========================================================================
    // SearchResult container
    // =========================================================================
    py::class_<crispritz::SearchResult>(m, "SearchResult",
                                        R"doc(
Aggregated result of a complete off-target search run.

Attributes
----------
off_targets : list[OffTarget]
    All deduplicated off-target hits across every queried chromosome.
profiles : list[PositionalProfile]
    One profile per input guide, in the same order as the input guide list.
)doc")
        .def_readonly("off_targets", &crispritz::SearchResult::off_targets)
        .def_readonly("profiles", &crispritz::SearchResult::profiles)
        .def("__len__", [](const crispritz::SearchResult& r) { return r.off_targets.size(); })
        .def("__repr__",
             [](const crispritz::SearchResult& r)
             {
                 return "<SearchResult hits=" + std::to_string(r.off_targets.size()) +
                        " guides=" + std::to_string(r.profiles.size()) + ">";
             });

    // =========================================================================
    // SearchParams struct
    // =========================================================================
    py::class_<crispritz::SearchParams>(m, "SearchParams",
                                        R"doc(
Parameters controlling a CRISPR off-target search run.

Parameters
----------
max_mismatches : int
    Maximum substitution mismatches allowed (>= 0).
max_dna_bulges : int
    Maximum DNA bulges allowed – gap in the guide / extra base in the target
    (>= 0).
max_rna_bulges : int
    Maximum RNA bulges allowed – gap in the target / extra base in the guide
    (>= 0).
pam_length : int
    Total PAM + guide pattern length (e.g. 23 for NGG with a 20-nt guide).
pam_limit : int
    PAM-only length (e.g. 3 for NGG).
pam_at_start : bool
    True for PAM-upstream systems (e.g. Cas12a ``TTTN``).
num_threads : int
    Number of OpenMP threads (default 1).
)doc")
        .def(py::init<>())
        .def_readwrite("max_mismatches", &crispritz::SearchParams::max_mismatches)
        .def_readwrite("max_dna_bulges", &crispritz::SearchParams::max_dna_bulges)
        .def_readwrite("max_rna_bulges", &crispritz::SearchParams::max_rna_bulges)
        .def_readwrite("pam_length", &crispritz::SearchParams::pam_length)
        .def_readwrite("pam_limit", &crispritz::SearchParams::pam_limit)
        .def_readwrite("pam_at_start", &crispritz::SearchParams::pam_at_start)
        .def_readwrite("num_threads", &crispritz::SearchParams::num_threads)
        .def("guide_length", &crispritz::SearchParams::guide_length,
             "Derived guide length: pam_length - pam_limit.");

    // =========================================================================
    // search_offtargets – primary entry point
    // =========================================================================
    m.def("search_offtargets",
          py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&,
                            const crispritz::SearchParams&>(&crispritz::search_offtargets),
          py::arg("bin_paths"), py::arg("guides"), py::arg("params"),
          R"doc(
Search pre-built TST index files for CRISPR off-target sites.

Parameters
----------
bin_paths : list[str]
    Paths to ``.bin`` partition files produced by ``build_tree``.
    All files for every chromosome to be searched must be provided.
guides : list[str]
    Guide RNA sequences WITHOUT the PAM, one per query.
    Example: ``["ACGTACGTACGTACGTACGT"]`` for a 20-nt guide.
params : SearchParams
    Search configuration (mismatches, bulges, PAM info, threading).

Returns
-------
SearchResult
    Contains ``off_targets`` (list of ``OffTarget``) and ``profiles``
    (list of ``PositionalProfile``, one per guide).

Raises
------
TSTSearchError
    If any parameter is invalid (negative limits, empty inputs, etc.).
TSTBuildError
    If a ``.bin`` file cannot be opened or has an unexpected format.
)doc");

    // Scalar-argument convenience overload.
    m.def("search_offtargets",
          py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&, int,
                            int, int, int, int, bool, int>(&crispritz::search_offtargets),
          py::arg("bin_paths"), py::arg("guides"), py::arg("max_mismatches"),
          py::arg("max_dna_bulges"), py::arg("max_rna_bulges"), py::arg("pam_length"),
          py::arg("pam_limit"), py::arg("pam_at_start"), py::arg("num_threads") = 1,
          R"doc(
Search pre-built TST index files for CRISPR off-target sites (scalar API).

This is a convenience overload that builds a ``SearchParams`` internally.

Parameters
----------
bin_paths : list[str]
    Paths to ``.bin`` partition files produced by ``build_tree``.
guides : list[str]
    Guide RNA sequences WITHOUT the PAM.
max_mismatches : int
    Maximum substitutions allowed (>= 0).
max_dna_bulges : int
    Maximum DNA bulges allowed (>= 0).
max_rna_bulges : int
    Maximum RNA bulges allowed (>= 0).
pam_length : int
    Total PAM + guide pattern length.
pam_limit : int
    PAM-only length.
pam_at_start : bool
    True for PAM-upstream systems (e.g. Cas12a).
num_threads : int, optional
    Number of OpenMP threads (default 1).

Returns
-------
SearchResult
    Contains ``off_targets`` and ``profiles``.

Raises
------
TSTSearchError
    If any parameter is invalid.
TSTBuildError
    If a ``.bin`` file cannot be opened or parsed.
)doc");
}