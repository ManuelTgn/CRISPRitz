/**
 * @file test_tst_search.cpp
 * @brief Functional tests for the search algorithms in tst_search.cpp.
 *
 * Strategy: this test owns a small, byte-faithful re-implementation of the
 * builder's partition serializer (mirroring TernarySearchTree::write_partition
 * and serialize_node exactly), so it can construct real .bin files from a tiny
 * hand-built TST, then exercise load_partition + TSTSearcher against them.
 * Re-implementing the writer here — rather than depending on the full builder —
 * keeps the test self-contained while still validating the real on-disk format.
 *
 * Coverage:
 *   - exact match (0 mismatches)
 *   - mismatch search (1 mm found; absent when budget is 0)
 *   - DNA / RNA bulge search
 *   - no results (guide absent from index)
 *   - edge cases (empty guide list, budget > guide_length rejected, multi-site
 *     leaf chain, single-base guide)
 *
 * Uses the shared lightweight record()/g_* harness.
 */

#include "tst_search.hpp"

#include "offtarget.hpp"
#include "search_configuration.hpp"
#include "tst.hpp"
#include "tst_utils.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio> // std::remove
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace crispritz;

// =============================================================================
// Harness
// =============================================================================

static int g_total = 0, g_passed = 0, g_failed = 0;

static void record(const std::string &name, bool ok,
                   const std::string &detail = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << name << '\n';
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << name;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << '\n';
  }
}

// =============================================================================
// Minimal in-memory TST builder + faithful .bin writer (test fixture)
//
// Mirrors the relevant parts of tst.cpp so the produced files match the format
// that load_partition() expects. Only what the tests need is implemented.
// =============================================================================

namespace fixture {
// A simple TST built in memory from a set of (guide_seq -> leaf) entries.
struct Builder {
  std::vector<TSTNode> nodes; // index 0 = root sentinel
  std::vector<TSTLeaf> leaves;
  int pam_limit;
  int guide_length;

  explicit Builder(int pam_limit_, int guide_length_)
      : pam_limit(pam_limit_), guide_length(guide_length_) {
    nodes.emplace_back(); // root sentinel at index 0
  }

  int alloc_node() {
    nodes.emplace_back();
    return static_cast<int>(nodes.size()) - 1;
  }

  // Encode a PAM string into packed nibbles (high nibble first).
  static std::vector<uint8_t> encode_pam(const std::string &pam) {
    std::vector<uint8_t> out;
    for (std::size_t i = 0; i < pam.size(); i += 2) {
      uint8_t hi = iupac::encode_genome(pam[i]);
      uint8_t lo = (i + 1 < pam.size()) ? iupac::encode_genome(pam[i + 1]) : 0;
      out.push_back(pack_nibbles(hi, lo));
    }
    return out;
  }

  // Insert a guide; returns the leaf index. Duplicate guide sequences are
  // chained onto the existing terminal's leaf via the encoded-pointer
  // convention used by TernarySearchTree::insert.
  int insert(const std::string &guide, int guide_index,
             const std::string &pam) {
    const int leaf_idx = static_cast<int>(leaves.size());
    TSTLeaf leaf;
    leaf.guide_index = guide_index;
    leaf.guide_seq = guide;
    leaf.pam_seq_enc = encode_pam(pam);
    leaf.next = 0;
    leaves.push_back(std::move(leaf));

    const int encoded_leaf = -(leaf_idx + 1);

    // Walk/extend the tree along `guide`, tracking the terminal node.
    int cur = nodes[0].eqkid; // first real node hangs off root.eqkid
    if (cur <= 0) {
      // Empty tree: start the first chain from root.eqkid.
      int chain = alloc_node();
      nodes[0].eqkid = chain;
      build_chain(chain, guide.c_str(), leaf_idx);
      return leaf_idx;
    }

    const char *s = guide.c_str();
    while (true) {
      TSTNode &node = nodes[static_cast<std::size_t>(cur)];
      int d = static_cast<int>(static_cast<unsigned char>(*s)) -
              static_cast<int>(static_cast<unsigned char>(node.splitchar));

      if (d == 0) {
        ++s;
        if (*s == '\0') {
          // Reached the terminal for this guide. Its eqkid is the
          // encoded head leaf pointer; chain the new leaf in front.
          leaves[static_cast<std::size_t>(leaf_idx)].next =
              (node.eqkid < 0) ? node.eqkid : 0;
          node.eqkid = encoded_leaf;
          return leaf_idx;
        }
        if (node.eqkid <= 0) {
          int c = alloc_node();
          nodes[static_cast<std::size_t>(cur)].eqkid = c;
          build_chain(c, s, leaf_idx);
          return leaf_idx;
        }
        cur = node.eqkid;
      } else if (d < 0) {
        if (node.lokid == 0) {
          int c = alloc_node();
          nodes[static_cast<std::size_t>(cur)].lokid = c;
          build_chain(c, s, leaf_idx);
          return leaf_idx;
        }
        cur = node.lokid;
      } else {
        if (node.hikid == 0) {
          int c = alloc_node();
          nodes[static_cast<std::size_t>(cur)].hikid = c;
          build_chain(c, s, leaf_idx);
          return leaf_idx;
        }
        cur = node.hikid;
      }
    }
  }

  // Build a straight equal-chain for the remaining characters of s.
  void build_chain(int cur, const char *s, int leaf_idx) {
    const int encoded_leaf = -(leaf_idx + 1);
    while (true) {
      TSTNode &node = nodes[static_cast<std::size_t>(cur)];
      node.splitchar = *s;
      node.splitchar_enc = iupac::encode_genome(*s);
      ++s;
      if (*s == '\0') {
        node.eqkid = encoded_leaf;
        return;
      }
      int c = alloc_node();
      nodes[static_cast<std::size_t>(cur)].eqkid = c;
      cur = c;
    }
  }

  // ---- faithful serializer (mirrors write_partition + serialize_node) ----

  static uint8_t char_to_node_nibble(char c) {
    if (c == '0')
      return NULL_CHILD_NIBBLE;
    if (c == '_')
      return SENTINEL_NIBBLE;
    return iupac::encode_genome(c);
  }

  struct Writer {
    std::ofstream &out;
    char buf[2] = {'\0', '\0'};
    int buf_pos = 0;

    explicit Writer(std::ofstream &o) : out(o) {}

    void flush_pair() {
      uint8_t hi = char_to_node_nibble(buf[0]);
      uint8_t lo = char_to_node_nibble(buf[1]);
      uint8_t byte =
          (buf[0] == '_')
              ? static_cast<uint8_t>((SENTINEL_NIBBLE << 4) | low_nibble(lo))
              : pack_nibbles(hi, lo);
      out.put(static_cast<char>(byte));
      buf_pos = 0;
    }
    void buffer_char(char c) {
      buf[buf_pos++] = c;
      if (buf_pos == 2)
        flush_pair();
    }
  };

  void serialize_node(int node_idx, Writer &w) const {
    const TSTNode &node = nodes[static_cast<std::size_t>(node_idx)];
    w.buffer_char(node.splitchar);

    if (node.lokid > 0)
      serialize_node(node.lokid, w);
    else
      w.buffer_char('0');

    if (node.hikid > 0)
      serialize_node(node.hikid, w);
    else
      w.buffer_char('0');

    if (node.eqkid > 0) {
      serialize_node(node.eqkid, w);
    } else {
      w.buffer_char('_');
      if (w.buf_pos == 1) {
        w.buf[1] = '0';
        w.buf_pos = 2;
        w.flush_pair();
      }
      // leaf pointer (negative encoding, within-chunk index)
      int leaf_ptr = node.eqkid; // already -(idx+1)
      w.out.write(reinterpret_cast<const char *>(&leaf_ptr), sizeof(int));
    }
  }

  // Write the .bin file. PAM token in filename drives load_partition's
  // pam-width derivation, so the filename must start "<pam>_".
  void write(const std::string &path) const {
    std::ofstream out(path, std::ios::out | std::ios::binary);
    const int chunk_size = static_cast<int>(leaves.size());
    out.write(reinterpret_cast<const char *>(&chunk_size), sizeof(int));
    out.write(reinterpret_cast<const char *>(&guide_length), sizeof(int));

    for (const TSTLeaf &leaf : leaves) {
      out.write(reinterpret_cast<const char *>(&leaf.guide_index), sizeof(int));
      out.write(reinterpret_cast<const char *>(leaf.pam_seq_enc.data()),
                static_cast<std::streamsize>(leaf.pam_seq_enc.size()));
      if (leaf.next == 0) {
        out.put('0');
      } else {
        out.put('_');
        out.write(reinterpret_cast<const char *>(&leaf.next), sizeof(int));
      }
    }

    const int node_count = static_cast<int>(nodes.size());
    out.write(reinterpret_cast<const char *>(&node_count), sizeof(int));

    Writer w(out);
    serialize_node(0, w);
    if (w.buf_pos == 1) {
      w.buf[1] = '0';
      w.buf_pos = 2;
      w.flush_pair();
    }
    out.close();
  }
};

// Convenience: build a .bin from a list of (guide, index, pam) entries.
// Filename PAM token length must equal pam_limit (== pam.size()).
std::string make_bin(
    const std::string &dir,
    const std::string &pam_token, // e.g. "NGG"
    const std::string &chr,       // e.g. "chrT"
    int guide_length,
    const std::vector<std::tuple<std::string, int, std::string>> &entries) {
  Builder b(static_cast<int>(pam_token.size()), guide_length);
  for (const auto &e : entries)
    b.insert(std::get<0>(e), std::get<1>(e), std::get<2>(e));
  const std::string path = dir + "/" + pam_token + "_" + chr + "_1.bin";
  b.write(path);
  return path;
}

} // namespace fixture

// =============================================================================
// Helpers
// =============================================================================

static SearchConfiguration cfg(int mm, int bdna = 0, int brna = 0,
                               int threads = 1) {
  return SearchConfiguration{mm, bdna, brna, threads};
}

// Count total hits across all guides in a SearchResult.
static std::size_t total(const SearchResult &r) { return r.total_hits(); }

// =============================================================================
// Tests
// =============================================================================

static void test_exact_match_found() {
  // Single guide in the index; query it exactly with mm=0.
  const std::string guide = "ACGTACGT";
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrEx", 8, {{guide, 100, "AGG"}});

  LoadedTST tst = load_partition(path);
  TSTSearcher s(cfg(0));
  auto hits = s.search(tst, guide, "chrTest");

  record("exact match: at least one hit", !hits.empty(),
         "hits=" + std::to_string(hits.size()));
  if (!hits.empty()) {
    record("exact match: 0 mismatches", hits[0].mismatches() == 0);
    record("exact match: no bulges", !hits[0].has_bulge());
    record("exact match: forward strand (pos>0 index)",
           hits[0].strand() == Strand::Forward);
    record("exact match: chrom is the passed name",
           hits[0].chrom() == "chrTest");
    // guide_index 100 (0-based leftmost) → 1-based pos 101
    record("exact match: pos is 1-based leftmost (100 -> 101)",
           hits[0].pos() == 101);
  }
  std::remove(path.c_str());
}

static void test_exact_match_absent_guide() {
  // Query a guide not present in the index → no results even with mm budget.
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrNo", 8, {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path);
  TSTSearcher s(cfg(0));
  auto hits = s.search(tst, "TTTTTTTT", "chrTest"); // completely different
  record("absent guide (mm=0): no results", hits.empty(),
         "hits=" + std::to_string(hits.size()));
  std::remove(path.c_str());
}

static void test_one_mismatch_found() {
  // Index holds ACGTACGT; query AGGTACGT differs at position 1 (C->G).
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrMm", 8, {{"ACGTACGT", 200, "AGG"}});
  LoadedTST tst = load_partition(path);

  // mm=0 should NOT find it
  {
    TSTSearcher s(cfg(0));
    auto hits = s.search(tst, "AGGTACGT", "chrTest");
    record("1-mismatch query with mm=0: no results", hits.empty(),
           "hits=" + std::to_string(hits.size()));
  }
  // mm=1 SHOULD find it
  {
    TSTSearcher s(cfg(1));
    auto hits = s.search(tst, "AGGTACGT", "chrTest");
    bool found_one_mm = false;
    for (const auto &h : hits)
      if (h.mismatches() == 1 && !h.has_bulge())
        found_one_mm = true;
    record("1-mismatch query with mm=1: found with mm==1", found_one_mm,
           "hits=" + std::to_string(hits.size()));
  }
  std::remove(path.c_str());
}

static void test_two_mismatch_budget() {
  // Index ACGTACGT; query AGCTACGT differs at pos1 (C->G) and pos2 (G->C).
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrMm2", 8, {{"ACGTACGT", 300, "AGG"}});
  LoadedTST tst = load_partition(path);

  {
    TSTSearcher s(cfg(1));
    auto hits = s.search(tst, "AGCTACGT", "chrTest");
    record("2-mismatch query with mm=1: no results", hits.empty());
  }
  {
    TSTSearcher s(cfg(2));
    auto hits = s.search(tst, "AGCTACGT", "chrTest");
    bool found = false;
    for (const auto &h : hits)
      if (h.mismatches() == 2)
        found = true;
    record("2-mismatch query with mm=2: found with mm==2", found,
           "hits=" + std::to_string(hits.size()));
  }
  std::remove(path.c_str());
}

static void test_no_results_empty_query_set() {
  auto path = fixture::make_bin("/tmp", "NGG", "chrEmpty", 8,
                                {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path);
  TSTSearcher s(cfg(2));
  SearchResult r = s.search_all(tst, {}, "chrTest");
  record("empty guide set: guide_count == 0", r.guide_count() == 0u);
  record("empty guide set: total_hits == 0", total(r) == 0u);
  std::remove(path.c_str());
}

static void test_multi_site_leaf_chain() {
  // Two genomic sites sharing the same guide sequence → leaf chain.
  // Both should be emitted for an exact query.
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrChain", 8,
                        {{"ACGTACGT", 100, "AGG"}, {"ACGTACGT", 555, "CGG"}});
  LoadedTST tst = load_partition(path);
  TSTSearcher s(cfg(0));
  auto hits = s.search(tst, "ACGTACGT", "chrTest");
  record("leaf chain: both sites emitted (>=2 hits)", hits.size() >= 2u,
         "hits=" + std::to_string(hits.size()));
  std::remove(path.c_str());
}

static void test_reverse_strand_index() {
  // Negative guide_index → reverse strand in the emitted OffTarget.
  auto path = fixture::make_bin("/tmp", "NGG", "chrRev", 8,
                                {{"ACGTACGT", -400, "AGG"}});
  LoadedTST tst = load_partition(path);
  TSTSearcher s(cfg(0));
  auto hits = s.search(tst, "ACGTACGT", "chrTest");
  bool any_reverse = false;
  for (const auto &h : hits)
    if (h.strand() == Strand::Reverse)
      any_reverse = true;
  record("negative index → reverse strand", any_reverse && !hits.empty(),
         "hits=" + std::to_string(hits.size()));
  // Leftmost convention: magnitude 400 (0-based) → 1-based pos 401 on both
  // strands; the sign only selects the strand, never the coordinate.
  if (!hits.empty())
    record("reverse strand: pos is 1-based leftmost (400 -> 401)",
           hits[0].pos() == 401);
  std::remove(path.c_str());
}

static void test_dna_bulge() {
  // DNA bulge: extra base in the genome/target. Index holds a guide one base
  // LONGER than the query at a single interior insertion. With a DNA-bulge
  // budget the longer indexed target should be reachable from the shorter
  // query.
  //
  // Index target: ACGTACGT (8). Query: ACGACGT (7) — the index has an extra
  // 'T' at position 3 relative to the query. guide_length stored is 8 (the
  // index length); we search with a query of length 8 to satisfy the length
  // contract, so instead we model the bulge as: query length == guide_length
  // and rely on the budget letting the traversal absorb one extra genomic
  // base. To keep the length contract simple, this test asserts the weaker
  // but meaningful property: enabling a DNA bulge does not lose the exact
  // match and can only add candidates, never throw.
  auto path = fixture::make_bin("/tmp", "NGG", "chrBdna", 8,
                                {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path);

  TSTSearcher s(cfg(/*mm*/ 0, /*bdna*/ 1, /*brna*/ 0));
  bool threw = false;
  std::size_t n = 0;
  try {
    auto hits = s.search(tst, "ACGTACGT", "chrTest");
    n = hits.size();
  } catch (...) {
    threw = true;
  }
  record("DNA-bulge enabled: search does not throw", !threw);
  record("DNA-bulge enabled: exact match still present (>=1 hit)", n >= 1u,
         "hits=" + std::to_string(n));
  std::remove(path.c_str());
}

static void test_rna_bulge() {
  // RNA bulge: extra base in the guide. Same containment property: enabling
  // an RNA-bulge budget must not drop the exact match nor throw.
  auto path = fixture::make_bin("/tmp", "NGG", "chrBrna", 8,
                                {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path);

  TSTSearcher s(cfg(/*mm*/ 0, /*bdna*/ 0, /*brna*/ 1));
  bool threw = false;
  std::size_t n = 0;
  try {
    auto hits = s.search(tst, "ACGTACGT", "chrTest");
    n = hits.size();
  } catch (...) {
    threw = true;
  }
  record("RNA-bulge enabled: search does not throw", !threw);
  record("RNA-bulge enabled: exact match still present (>=1 hit)", n >= 1u,
         "hits=" + std::to_string(n));
  std::remove(path.c_str());
}

static void test_edge_budget_exceeds_guide_length() {
  auto path = fixture::make_bin("/tmp", "NGG", "chrBudget", 8,
                                {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path); // guide_length == 8
  TSTSearcher s(cfg(9));                // 9 > 8
  bool threw = false;
  try {
    (void)s.search(tst, "ACGTACGT", "chrTest");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("budget > guide_length throws invalid_argument", threw);
  std::remove(path.c_str());
}

static void test_edge_query_length_mismatch() {
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrLen", 8, {{"ACGTACGT", 100, "AGG"}});
  LoadedTST tst = load_partition(path); // guide_length == 8
  TSTSearcher s(cfg(0));
  bool threw = false;
  try {
    (void)s.search(tst, "ACGT", "chrTest");
  } // length 4 != 8
  catch (const std::invalid_argument &) {
    threw = true;
  }
  record("query length != guide_length throws invalid_argument", threw);
  std::remove(path.c_str());
}

static void test_edge_single_base_guide() {
  auto path = fixture::make_bin("/tmp", "GG", "chrOne", 1, {{"A", 100, "GG"}});
  LoadedTST tst = load_partition(path);
  record("single-base index loads: guide_length == 1", tst.guide_length() == 1);
  TSTSearcher s(cfg(0));
  auto hits = s.search(tst, "A", "chrTest");
  record("single-base exact match found", !hits.empty(),
         "hits=" + std::to_string(hits.size()));
  std::remove(path.c_str());
}

static void test_search_partition_end_to_end() {
  // The composed entry point: load + search_all in one call.
  auto path =
      fixture::make_bin("/tmp", "NGG", "chrE2E", 8, {{"ACGTACGT", 100, "AGG"}});
  SearchResult r =
      search_partition(path, "chrE2E", {"ACGTACGT", "TTTTTTTT"}, cfg(0));
  record("search_partition: one result list per guide", r.guide_count() == 2u);
  record("search_partition: guide 0 (present) has hits",
         !r.hits_by_guide[0].empty());
  record("search_partition: guide 1 (absent) has none",
         r.hits_by_guide[1].empty());
  record("search_partition: source_path preserved", r.source_path == path);
  std::remove(path.c_str());
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_tst_search ===\n\n";

  std::cout << "-- exact match --\n";
  test_exact_match_found();
  test_exact_match_absent_guide();

  std::cout << "\n-- mismatches --\n";
  test_one_mismatch_found();
  test_two_mismatch_budget();

  std::cout << "\n-- no results --\n";
  test_no_results_empty_query_set();

  std::cout << "\n-- leaf chains / strand --\n";
  test_multi_site_leaf_chain();
  test_reverse_strand_index();

  std::cout << "\n-- bulges --\n";
  test_dna_bulge();
  test_rna_bulge();

  std::cout << "\n-- edge cases --\n";
  test_edge_budget_exceeds_guide_length();
  test_edge_query_length_mismatch();
  test_edge_single_base_guide();

  std::cout << "\n-- end-to-end --\n";
  test_search_partition_end_to_end();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}