/**
 * @file test_pam_search.cpp
 * @brief Unit tests for pam_search.cpp
 *
 * Covers:
 *  - CompactGenome construction and random access
 *  - search_pam_sites / search_pam_sites_fast (PAM-at-end, PAM-at-start)
 *  - Exact matches on forward and reverse strands
 *  - IUPAC ambiguity codes in the PAM pattern
 *  - Out-of-bound guard (sites too close to the chromosome edges)
 *  - Empty genome and empty result edge cases
 *  - Multiple-thread consistency (results with 1 vs N threads must be equal)
 */

#include "nucleotide_encoding.hpp"
#include "pam_search.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using pam::CompactGenome;
using pam::NucleotideEncoder;
using pam::search_pam_sites;
using pam::search_pam_sites_fast;
using pam::SearchParams;

// -----------------------------------------------------------------------------
// Minimal test harness (same convention as test_nucleotide_encoding.cpp)
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
// Helpers
// -----------------------------------------------------------------------------

/** @brief Return a sorted copy of a vector<int>. */
static std::vector<int> sorted(std::vector<int> v) {
  std::sort(v.begin(), v.end());
  return v;
}

/** @brief True when every element of 'expected' appears in 'actual'. */
static bool contains_all(const std::vector<int> &actual,
                         const std::vector<int> &expected) {
  std::set<int> s(actual.begin(), actual.end());
  for (int x : expected)
    if (s.find(x) == s.end())
      return false;
  return true;
}

// -----------------------------------------------------------------------------
// CompactGenome tests
// -----------------------------------------------------------------------------

/** @brief Constructor must store the sequence length correctly. */
static void test_compact_genome_size() {
  const std::string seq = "ACGTACGT";
  CompactGenome cg(seq);
  record("CompactGenome size == seq length", cg.size() == seq.size(),
         "got " + std::to_string(cg.size()));
}

/** @brief Operator[] must return the same encoding as
 * NucleotideEncoder::encode_genome. */
static void test_compact_genome_access() {
  const std::string seq = "ACGTNRYMKSWHBVD";
  CompactGenome cg(seq);
  bool ok = true;
  for (size_t i = 0; i < seq.size(); ++i) {
    uint8_t expected = NucleotideEncoder::encode_genome(seq[i]);
    uint8_t got = cg[i];
    if (got != expected) {
      ok = false;
      std::cout << "    mismatch at " << i << ": char=" << seq[i]
                << " expected=" << (int)expected << " got=" << (int)got << "\n";
    }
  }
  record("CompactGenome operator[] matches encode_genome", ok);
}

/** @brief data() pointer must be non-null and bytes() must be ceil(size/2). */
static void test_compact_genome_data_pointer() {
  const std::string seq = "ACGTACGT"; // 8 chars -> 4 bytes
  CompactGenome cg(seq);
  record("CompactGenome data() != nullptr", cg.data() != nullptr);
  record("CompactGenome bytes() == ceil(size/2)",
         cg.bytes() == (seq.size() + 1) / 2,
         "bytes=" + std::to_string(cg.bytes()));
}

/** @brief Size 1 genome (edge of even/odd packing). */
static void test_compact_genome_single_base() {
  CompactGenome cg("A");
  record("CompactGenome single base size == 1", cg.size() == 1u);
  record("CompactGenome single base [0] == encode_genome(A)",
         cg[0] == NucleotideEncoder::encode_genome('A'));
}

// -----------------------------------------------------------------------------
// search_pam_sites – PAM at end (pam_at_start = false)
// -----------------------------------------------------------------------------

/**
 * @brief Hand-crafted genome with a known NGG at position 5 on the forward
 *        strand.  With pam_length=3, pam_limit=3, the PAM itself occupies
 *        positions [5,7].  The guide start position returned must be 0
 *        (guide occupies [0,4] but with pam_length=5 that would be 0).
 *
 *        Here we keep pam_length == pam_limit == 3 (guide_length = 0) as a
 *        minimal smoke-test focused purely on the position arithmetic.
 *
 *        Genome (20 bases): AAAAATGGAAAAAAAAAAAA
 *                            0123456789...
 *        NGG at position 5 (0-based).  pam_length=3, pam_limit=3.
 *        Expected positive hit (guide_start): 5 + 3 - 1 - (3 - 1) = 5
 */
static void test_pam_at_end_forward_hit() {
  // Genome: NGG at index 5
  const std::string genome = "AAAAATGGAAAAAAAAAAAA"; // 20 chars
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites = search_pam_sites("NGG", genome, params);
  // At least one positive (forward) site must be found
  bool any_positive = false;
  for (int s : sites)
    if (s >= 0) {
      any_positive = true;
      break;
    }
  record("pam_at_end: forward hit found", any_positive,
         "sites=" + std::to_string(sites.size()));
}

/**
 * @brief Reverse strand: place CCN (rc of NGG) so a reverse-strand hit is
 *        detected.
 *
 *        Genome: AAAAACCGAAAAAAAAAAAA
 *        CCN at index 5.  RC of CCN is NGG, so a reverse hit is expected.
 *        Negative site indices indicate reverse-strand matches.
 */
static void test_pam_at_end_reverse_hit() {
  const std::string genome = "AAAAACCGAAAAAAAAAAAA"; // 20 chars
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites = search_pam_sites("NGG", genome, params);
  bool any_negative = false;
  for (int s : sites)
    if (s < 0) {
      any_negative = true;
      break;
    }
  record("pam_at_end: reverse hit found", any_negative,
         "sites=" + std::to_string(sites.size()));
}

/**
 * @brief No PAM in the sequence -> result must be empty.
 */
static void test_pam_at_end_no_hit() {
  const std::string genome = "AAAAAAAAAAAAAAAAAAAA"; // 20 A's, no GG
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites = search_pam_sites("NGG", genome, params);
  record("pam_at_end: no hit returns empty", sites.empty(),
         "sites=" + std::to_string(sites.size()));
}

// -----------------------------------------------------------------------------
// search_pam_sites – PAM at start (pam_at_start = true)
// -----------------------------------------------------------------------------

/**
 * @brief PAM-at-start (Cas12a-style): TTT at position 3 should register a
 *        forward (negative index by convention) hit.
 *
 *        Genome: AAATTTAAAAAAAAAAAAAAA (21 chars)
 *        TTT at index 3.  pam_length=20, pam_limit=3, pam_at_start=true.
 */
static void test_pam_at_start_forward_hit() {
  const std::string genome = "AAATTTAAAAAAAAAAAAAAA"; // 21 chars
  // pam_length=6 (3 PAM + 3 placeholder guide), pam_limit=3
  SearchParams params(6, 3, /*pam_at_start=*/true, 1);

  auto sites = search_pam_sites("TTT", genome, params);
  record("pam_at_start: at least one hit found", !sites.empty(),
         "sites=" + std::to_string(sites.size()));
}

// -----------------------------------------------------------------------------
// search_pam_sites_fast – consistency with search_pam_sites
// -----------------------------------------------------------------------------

/**
 * @brief search_pam_sites_fast must return the same set of sites as
 *        search_pam_sites for the same genome and PAM.
 */
static void test_fast_vs_normal_consistency() {
  const std::string genome =
      "ACGTACGTNGGACGTACGTNGGACGTACGT" // 30 chars, NGG at 8 and 18
      "ACGTACGTACGTACGTACGT";          // + 20 = 50 total
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites_normal = search_pam_sites("NGG", genome, params);
  CompactGenome cg(genome);
  auto sites_fast = search_pam_sites_fast("NGG", cg, params);

  record("fast vs normal: same count", sites_normal.size() == sites_fast.size(),
         "normal=" + std::to_string(sites_normal.size()) +
             " fast=" + std::to_string(sites_fast.size()));

  auto sn = sorted(sites_normal);
  auto sf = sorted(sites_fast);
  record("fast vs normal: same elements", sn == sf);
}

// -----------------------------------------------------------------------------
// IUPAC ambiguity in the PAM pattern
// -----------------------------------------------------------------------------

/**
 * @brief PAM 'NRG' (R = A|G) must match 'NAG' and 'NGG' on the genome.
 *
 *        Genome: AAAAATAGTGGAAAAAAAAAA (21 chars)
 *                           ^  ^
 *                      NAG@5   NGG@8
 */
static void test_iupac_pam_ambiguity() {
  const std::string genome = "AAAAATAGTGGAAAAAAAAAA"; // 21 chars
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites = search_pam_sites("NRG", genome, params);
  // We expect at least two forward hits (NAG and NGG satisfy NRG)
  long pos_count =
      std::count_if(sites.begin(), sites.end(), [](int s) { return s >= 0; });
  record("IUPAC NRG matches NAG and NGG: >= 2 fwd hits", pos_count >= 2,
         "positive hits=" + std::to_string(pos_count));
}

/**
 * @brief PAM 'NYG' (Y = C|T) must NOT match 'NGG' on the genome.
 *
 *        Genome: AAAAAGGGAAAAAAAAAAAAA (21 chars)
 *        NYG requires the second position to be C or T; G fails that check.
 */
static void test_iupac_pam_no_match() {
  const std::string genome = "AAAAAGGGAAAAAAAAAAAAA"; // 21 chars
  SearchParams params(3, 3, /*pam_at_start=*/false, 1);

  auto sites_nyg = search_pam_sites("NYG", genome, params);
  // No forward hits expected (NGG does not match NYG)
  long pos_count = std::count_if(sites_nyg.begin(), sites_nyg.end(),
                                 [](int s) { return s >= 0; });
  record("IUPAC NYG does not match NGG on fwd strand", pos_count == 0,
         "positive hits=" + std::to_string(pos_count));
}

// -----------------------------------------------------------------------------
// Multi-thread consistency
// -----------------------------------------------------------------------------

/**
 * @brief Running search_pam_sites_fast with 1 thread and with 4 threads must
 *        produce the same set of results.
 */
static void test_multithread_consistency() {
  // Build a genome long enough that OMP chunks differ between 1 and 4 threads
  std::string genome;
  genome.reserve(10000);
  const char bases[] = "ACGT";
  for (int i = 0; i < 10000; ++i)
    genome += bases[i % 4];

  // Inject a handful of known NGG sites
  for (int pos : {100, 1004, 3000, 7500, 9990}) {
    if (pos + 3 <= (int)genome.size()) {
      genome[pos] = 'A'; // N placeholder
      genome[pos + 1] = 'G';
      genome[pos + 2] = 'G';
    }
  }

  CompactGenome cg(genome);
  SearchParams p1(3, 3, false, 1);
  SearchParams p4(3, 3, false, 4);

  auto s1 = sorted(search_pam_sites_fast("NGG", cg, p1));
  auto s4 = sorted(search_pam_sites_fast("NGG", cg, p4));

  record("multithread: 1-thread count == 4-thread count",
         s1.size() == s4.size(),
         "1thr=" + std::to_string(s1.size()) +
             " 4thr=" + std::to_string(s4.size()));
  record("multithread: 1-thread sites == 4-thread sites", s1 == s4);
}

// -----------------------------------------------------------------------------
// Edge cases
// -----------------------------------------------------------------------------

/**
 * @brief A genome that is exactly pam_limit characters long.  There is no room
 *        for a guide, so positions must be filtered by the bound check.  We
 *        expect either zero results or only results satisfying the guard.
 */
static void test_minimal_genome() {
  const std::string genome = "NGG"; // exactly pam_limit=3 chars
  SearchParams params(3, 3, false, 1);
  // Should not crash; result is empty or one site (guide_start = 0)
  bool threw = false;
  std::vector<int> sites;
  try {
    sites = search_pam_sites("NGG", genome, params);
  } catch (...) {
    threw = true;
  }
  record("minimal genome does not throw", !threw);
  // Any returned positive site must be >= 0
  bool valid = true;
  for (int s : sites)
    if (s < 0) {
      valid = false;
      break;
    }
  record("minimal genome: no invalid (negative) forward sites", valid);
}

/**
 * @brief Genome made of all N's: the encoded value is 0b0000, which ANDs to 0
 *        with any PAM code.  No PAM hit should be found.
 */
static void test_all_n_genome() {
  const std::string genome(50, 'N');
  SearchParams params(3, 3, false, 1);
  auto sites = search_pam_sites("NGG", genome, params);
  // N (0b0000) & G (0b0100) == 0 -> no match
  record("all-N genome: no hits", sites.empty(),
         "sites=" + std::to_string(sites.size()));
}

/**
 * @brief Searching with 0 threads must not crash (constructor clamps to 1).
 */
static void test_zero_threads_no_crash() {
  const std::string genome = "AAAAANGGAAAAAAAAAAAA";
  SearchParams params(3, 3, false, 0); // 0 -> clamped to 1
  bool threw = false;
  try {
    search_pam_sites("NGG", genome, params);
  } catch (...) {
    threw = true;
  }
  record("0 threads does not throw", !threw);
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------

int main() {
  std::cout << "=== test_pam_search ===\n\n";

  std::cout << "-- CompactGenome --\n";
  test_compact_genome_size();
  test_compact_genome_access();
  test_compact_genome_data_pointer();
  test_compact_genome_single_base();

  std::cout << "\n-- PAM at end (forward) --\n";
  test_pam_at_end_forward_hit();
  test_pam_at_end_reverse_hit();
  test_pam_at_end_no_hit();

  std::cout << "\n-- PAM at start --\n";
  test_pam_at_start_forward_hit();

  std::cout << "\n-- fast vs normal consistency --\n";
  test_fast_vs_normal_consistency();

  std::cout << "\n-- IUPAC ambiguity --\n";
  test_iupac_pam_ambiguity();
  test_iupac_pam_no_match();

  std::cout << "\n-- multi-thread consistency --\n";
  test_multithread_consistency();

  std::cout << "\n-- edge cases --\n";
  test_minimal_genome();
  test_all_n_genome();
  test_zero_threads_no_crash();

  std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}
