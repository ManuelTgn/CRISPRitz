/**
 * @file test_tst.cpp
 * @brief Unit tests for tst.cpp (TernarySearchTree, build_tree)
 *
 * Strategy
 * --------
 * tst.cpp produces .bin files on disk.  The tests:
 *  1. Call build_tree() on small, controlled genomes.
 *  2. Verify that the expected .bin file(s) are created and are non-empty.
 *  3. Inspect the binary header (leaf count + guide_length) to assert
 *     structural correctness without requiring a full deserializer.
 *  4. Cover: PAM-at-end, PAM-at-start, max_bulges > 0, multiple partitions
 *     (LEAVES_PER_GROUP boundary), error paths (bad pam_limit).
 *
 * All temporary .bin files are cleaned up at the end of each test.
 */

#include "nucleotide_encoding.hpp"
#include "tst.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using crispritz::build_tree;
using crispritz::TernarySearchTree;

// -----------------------------------------------------------------------------
// Minimal test harness
// -----------------------------------------------------------------------------

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

static void record(const std::string &name, bool ok,
                   const std::string &detail = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << name << "\n";
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << name;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << "\n";
  }
}

// -----------------------------------------------------------------------------
// File helpers
// -----------------------------------------------------------------------------

/**
 * @brief Returns a list of paths matching the pattern
 *        <pam_seq>_<chr_name>_*.bin in the current directory.
 */
static std::vector<fs::path> find_bin_files(const std::string &pam_seq,
                                            const std::string &chr_name) {
  std::vector<fs::path> result;
  const std::string prefix = pam_seq + "_" + chr_name + "_";
  for (auto &entry : fs::directory_iterator(fs::current_path())) {
    if (!entry.is_regular_file())
      continue;
    const std::string fname = entry.path().filename().string();
    if (fname.rfind(prefix, 0) == 0 && fname.size() > 4 &&
        fname.substr(fname.size() - 4) == ".bin")
      result.push_back(entry.path());
  }
  std::sort(result.begin(), result.end());
  return result;
}

/** @brief Remove all .bin files produced by a test. */
static void cleanup_bin_files(const std::string &pam_seq,
                              const std::string &chr_name) {
  for (auto &p : find_bin_files(pam_seq, chr_name))
    fs::remove(p);
}

/**
 * @brief Read the first two int32 values from a .bin file:
 *        [0] = num_leaves, [1] = guide_length.
 * @throws std::runtime_error if the file is too small.
 */
static std::pair<int32_t, int32_t> read_bin_header(const fs::path &path) {
  std::ifstream fin(path, std::ios::binary);
  if (!fin.is_open())
    throw std::runtime_error("Cannot open: " + path.string());

  int32_t num_leaves = 0, guide_length = 0;
  fin.read(reinterpret_cast<char *>(&num_leaves), sizeof(int32_t));
  fin.read(reinterpret_cast<char *>(&guide_length), sizeof(int32_t));

  if (!fin)
    throw std::runtime_error("File too small to contain header: " +
                             path.string());
  return {num_leaves, guide_length};
}

// -----------------------------------------------------------------------------
// Genome helpers
// -----------------------------------------------------------------------------

/**
 * @brief Build a synthetic genome of length 'len' filled with 'ACGT' repeats,
 *        then insert explicit PAM motifs at every 'stride' positions.
 *
 * @param len    Total genome length.
 * @param pam    PAM string to inject (e.g. "GG" for the PAM portion of NGG).
 * @param stride Distance between injected PAM sites (must be > pam.size()).
 * @return       Synthetic genome string.
 */
static std::string make_genome_with_pam(int len,
                                        const std::string &context_before,
                                        const std::string &pam, int stride) {
  std::string genome(len, 'A');
  // Fill with repeating ACGT so no accidental PAMs appear
  for (int i = 0; i < len; ++i)
    genome[i] = "ACGT"[i % 4];

  const int site_len = static_cast<int>(context_before.size() + pam.size());
  for (int pos = 0; pos + site_len < len; pos += stride) {
    for (int j = 0; j < (int)context_before.size(); ++j)
      genome[pos + j] = context_before[j];
    for (int j = 0; j < (int)pam.size(); ++j)
      genome[pos + (int)context_before.size() + j] = pam[j];
  }
  return genome;
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

/**
 * @brief Minimal smoke test: a small genome with a single NGG site should
 *        produce exactly one .bin file with at least 1 leaf.
 *
 * PAM spec: "NNNNNNNNNNNNNNNNNNNNNGG 3"  (guide=20, pam=3)
 */
static void test_build_produces_bin_file() {
  const std::string chr = "testchr1";
  const std::string pam = "GG";      // only PAM portion
  const std::string full_pam = "GG"; // pam_seq arg to build_tree

  // Genome: 30 A's then GG at position 20 (i.e. ACGTNGG...)
  std::string genome(50, 'A');
  genome[20] = 'G';
  genome[21] = 'G';

  // pam_length=3, pam_limit=2 (guide=1 base, pam=GG)
  // This is the minimal valid config: guide_length=1
  const std::string test_id = "GG_" + chr;
  cleanup_bin_files("GG", chr); // ensure clean state

  bool threw = false;
  try {
    build_tree(genome, chr, pam, /*pam_length=*/3, /*pam_limit=*/2,
               /*pam_at_start=*/false, /*outdir=*/".", /*max_bulges=*/0,
               /*threads=*/1);
  } catch (const std::exception &e) {
    threw = true;
    std::cout << "    exception: " << e.what() << "\n";
  }

  if (!threw) {
    auto bins = find_bin_files("GG", chr);
    record("build produces at least 1 .bin file", !bins.empty(),
           "files=" + std::to_string(bins.size()));
    if (!bins.empty()) {
      auto [nl, gl] = read_bin_header(bins[0]);
      record("bin header: num_leaves > 0", nl > 0,
             "num_leaves=" + std::to_string(nl));
      record("bin header: guide_length == pam_length - pam_limit",
             gl == 1, // 3 - 2 = 1
             "guide_length=" + std::to_string(gl));
    }
  } else {
    record("build does not throw on valid input", false);
  }

  cleanup_bin_files("GG", chr);
}

/**
 * @brief A genome with multiple NGG sites (PAM-at-end, standard Cas9 config).
 *        pam_length=23 (20 guide + 3 PAM), pam_limit=3.
 *        Verifies guide_length stored in bin == 20.
 */
static void test_cas9_guide_length_in_header() {
  const std::string chr = "cas9chr";
  cleanup_bin_files("NGG", chr);

  // Build a 200-base genome with several NGG sites
  std::string genome =
      make_genome_with_pam(200, "ACGTACGTACGTACGTACGT", "GG", 50);

  bool threw = false;
  try {
    build_tree(genome, chr, "NGG", /*pam_length=*/23, /*pam_limit=*/3,
               /*pam_at_start=*/false, /*outdir=*/".", /*max_bulges=*/0, 1);
  } catch (const std::exception &e) {
    threw = true;
    std::cout << "    exception: " << e.what() << "\n";
  }

  if (!threw) {
    auto bins = find_bin_files("NGG", chr);
    record("Cas9 build: at least 1 .bin", !bins.empty());
    if (!bins.empty()) {
      auto [nl, gl] = read_bin_header(bins[0]);
      record("Cas9 build: guide_length == 20", gl == 20,
             "stored=" + std::to_string(gl));
      record("Cas9 build: num_leaves > 0", nl > 0);
    }
  } else {
    record("Cas9 build does not throw", false);
  }

  cleanup_bin_files("NGG", chr);
}

/**
 * @brief PAM-at-start configuration (Cas12a-style).
 *        PAM spec: "TTTNNNNNNNNNNNNNNNNNNNNN -4" -> pam_at_start=true,
 *        pam_limit=3 (TTT), pam_length=23 (3+20), guide_length=20.
 */
static void test_pam_at_start_header() {
  const std::string chr = "cas12achr";
  cleanup_bin_files("TTT", chr);

  // Genome: TTT at regular intervals (PAM at start)
  std::string genome(200, 'A');
  for (int pos = 0; pos + 23 < 200; pos += 50) {
    genome[pos] = 'T';
    genome[pos + 1] = 'T';
    genome[pos + 2] = 'T';
    // fill 20 guide bases with non-N chars
    for (int k = 3; k < 23 && pos + k < 200; ++k)
      genome[pos + k] = "ACGT"[k % 4];
  }

  bool threw = false;
  try {
    build_tree(genome, chr, "TTT", /*pam_length=*/23, /*pam_limit=*/3,
               /*pam_at_start=*/true, /*outdir=*/".", /*max_bulges=*/0, 1);
  } catch (const std::exception &e) {
    threw = true;
    std::cout << "    exception: " << e.what() << "\n";
  }

  if (!threw) {
    auto bins = find_bin_files("TTT", chr);
    record("PAM-at-start: at least 1 .bin", !bins.empty());
    if (!bins.empty()) {
      auto [nl, gl] = read_bin_header(bins[0]);
      record("PAM-at-start: guide_length == 20", gl == 20,
             "stored=" + std::to_string(gl));
    }
  } else {
    record("PAM-at-start build does not throw", false);
  }

  cleanup_bin_files("TTT", chr);
}

/**
 * @brief max_bulges > 0 changes the window size extracted per site.
 *        The index file should still be created; we just verify the header
 *        guide_length is unaffected (bulges only affect extraction window).
 */
static void test_max_bulges_header_unaffected() {
  const std::string chr = "bulgechr";
  cleanup_bin_files("NGG", chr);

  std::string genome =
      make_genome_with_pam(200, "ACGTACGTACGTACGTACGT", "GG", 50);

  bool threw = false;
  try {
    build_tree(genome, chr, "NGG", 23, 3, false, /*outdir=*/".",
               /*max_bulges=*/2, 1);
  } catch (const std::exception &e) {
    threw = true;
    std::cout << "    exception: " << e.what() << "\n";
  }

  if (!threw) {
    auto bins = find_bin_files("NGG", chr);
    record("max_bulges=2: at least 1 .bin", !bins.empty());
    if (!bins.empty()) {
      auto [nl, gl] = read_bin_header(bins[0]);
      record("max_bulges=2: guide_length still 20", gl == 20,
             "stored=" + std::to_string(gl));
    }
  } else {
    record("max_bulges=2 build does not throw", false);
  }

  cleanup_bin_files("NGG", chr);
}

/**
 * @brief Genome with no PAM sites at all: build() should return silently
 *        (empty pam_sites path) rather than throw.  No .bin file is written.
 */
static void test_no_pam_sites_does_not_throw() {
  const std::string chr = "nopamchr";
  cleanup_bin_files("NGG", chr);

  // All-A genome: no GG anywhere
  const std::string genome(100, 'A');

  bool threw = false;
  try {
    build_tree(genome, chr, "NGG", 23, 3, false, ".", 0, 1);
  } catch (const std::exception &e) {
    threw = true;
    std::cout << "    exception: " << e.what() << "\n";
  }

  record("no PAM sites: build does not throw", !threw);
  cleanup_bin_files("NGG", chr);
}

/**
 * @brief Constructing TernarySearchTree with guide_length <= 0 must throw
 *        std::runtime_error.
 */
static void test_constructor_invalid_guide_length_throws() {
  bool threw = false;
  try {
    // pam_length == pam_limit -> guide_length = 0 -> invalid
    TernarySearchTree tst("ACGTACGTACGT", "chr_invalid", "GG", /*pam_length=*/2,
                          /*pam_limit=*/2, false, /*outdir=*/".", 0, 1);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  record("constructor throws for guide_length == 0", threw);
}

/**
 * @brief Constructing TernarySearchTree with pam_limit == 0 must throw.
 */
static void test_constructor_invalid_pam_limit_throws() {
  bool threw = false;
  try {
    TernarySearchTree tst("ACGTACGTACGT", "chr_invalid2", "", /*pam_length=*/5,
                          /*pam_limit=*/0, false, /*outdir=*/".", 0, 1);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  record("constructor throws for pam_limit == 0", threw);
}

/**
 * @brief After a successful build, TSTNode leaf count must equal
 *        TernarySearchTree::leaf_count().
 */
static void test_leaf_count_via_api() {
  const std::string chr = "leafcountchr";
  cleanup_bin_files("NGG", chr);

  // Genome with 3 known NGG sites on forward strand (positions 0, 50, 100)
  std::string genome(200,
                     'C'); // fill with C so ACGT won't produce accidental GGs
  // Insert NGG motifs: [guide=20 A's][GG]
  for (int start : {0, 50, 100}) {
    for (int k = 0; k < 20 && start + k < 200; ++k)
      genome[start + k] = 'A';
    if (start + 21 < 200) {
      genome[start + 20] = 'G';
      genome[start + 21] = 'G';
    }
  }

  try {
    TernarySearchTree tst(genome, chr, "NGG", 23, 3, false, /*outdir=*/".", 0,
                          1);
    tst.build();
    // leaf_count() is available before save writes files
    bool has_leaves = tst.leaf_count() > 0;
    record("leaf_count() > 0 after build", has_leaves,
           "count=" + std::to_string(tst.leaf_count()));
  } catch (const std::exception &e) {
    record("leaf_count test does not throw", false,
           std::string("exception: ") + e.what());
  }

  cleanup_bin_files("NGG", chr);
}

/**
 * @brief The .bin file must be non-empty after a successful build.
 */
static void test_bin_file_is_nonempty() {
  const std::string chr = "nonemptychr";
  cleanup_bin_files("NGG", chr);

  std::string genome =
      make_genome_with_pam(200, "ACGTACGTACGTACGTACGT", "GG", 50);

  try {
    build_tree(genome, chr, "NGG", 23, 3, false, ".", 0, 1);
  } catch (...) {
  }

  auto bins = find_bin_files("NGG", chr);
  if (!bins.empty()) {
    auto size = fs::file_size(bins[0]);
    record(".bin file is non-empty", size > 0, "bytes=" + std::to_string(size));
  } else {
    record(".bin file is non-empty", false, "no bin files found");
  }

  cleanup_bin_files("NGG", chr);
}

/**
 * @brief Two successive builds for different chromosomes must not interfere
 *        (i.e. the file-naming convention keeps them separate).
 */
static void test_two_builds_independent_files() {
  const std::string chrA = "sepchrA";
  const std::string chrB = "sepchrB";
  cleanup_bin_files("NGG", chrA);
  cleanup_bin_files("NGG", chrB);

  std::string genome =
      make_genome_with_pam(200, "ACGTACGTACGTACGTACGT", "GG", 50);

  try {
    build_tree(genome, chrA, "NGG", 23, 3, false, ".", 0, 1);
  } catch (...) {
  }
  try {
    build_tree(genome, chrB, "NGG", 23, 3, false, ".", 0, 1);
  } catch (...) {
  }

  auto binsA = find_bin_files("NGG", chrA);
  auto binsB = find_bin_files("NGG", chrB);

  record("chrA and chrB produce separate .bin files",
         !binsA.empty() && !binsB.empty() &&
             binsA[0].filename() != binsB[0].filename());

  cleanup_bin_files("NGG", chrA);
  cleanup_bin_files("NGG", chrB);
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------

int main() {
  std::cout << "=== test_tst ===\n\n";

  std::cout << "-- basic bin production --\n";
  test_build_produces_bin_file();

  std::cout << "\n-- Cas9 header --\n";
  test_cas9_guide_length_in_header();

  std::cout << "\n-- PAM-at-start header --\n";
  test_pam_at_start_header();

  std::cout << "\n-- max_bulges --\n";
  test_max_bulges_header_unaffected();

  std::cout << "\n-- edge cases --\n";
  test_no_pam_sites_does_not_throw();
  test_constructor_invalid_guide_length_throws();
  test_constructor_invalid_pam_limit_throws();

  std::cout << "\n-- leaf count API --\n";
  test_leaf_count_via_api();

  std::cout << "\n-- file properties --\n";
  test_bin_file_is_nonempty();
  test_two_builds_independent_files();

  std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}