/**
 * @file test_offtarget.cpp
 * @brief Unit tests for the OffTarget domain model (offtarget.hpp /
 * offtarget.cpp).
 *
 * Test surface:
 *   - strand_from_char()          — round-trip and rejection of bad input
 *   - to_char()                   — constexpr conversion
 *   - OffTarget construction      — valid paths and every validation failure
 *   - Accessors                   — values match constructor arguments
 *   - total_edit_distance()       — additive formula
 *   - has_bulge()                 — derived from bulge counts
 *   - bulge_type()                — all four classification strings
 *   - operator==  / operator!=    — locus-based identity (excludes edit counts)
 *   - operator<                   — chrom → pos → strand → total_edit_distance
 *   - std::hash<OffTarget>        — consistent with operator==; usable in sets
 *   - to_tsv_row()                — column values, field count, custom
 * separator
 *   - tsv_header()                — column count and canonical names
 *   - locus_string()              — "chrom:pos(strand)" format
 *   - Edge cases                  — min position, long names, zero edit
 * distance
 *
 * Uses the same lightweight harness (record + g_* counters) as the existing
 * TST tests so it can be registered with the project's add_crispritz_test macro
 * without pulling in an external framework.
 */

#include "offtarget.hpp"

#include <algorithm> // std::sort
#include <iostream>
#include <limits>    // std::numeric_limits
#include <sstream>   // std::ostringstream (stream operator test)
#include <stdexcept> // std::invalid_argument
#include <string>
#include <unordered_set>
#include <vector>

using crispritz::OffTarget;
using crispritz::Strand;
using crispritz::strand_from_char;
using crispritz::to_char;

// =============================================================================
// Minimal test harness — mirrors existing CRISPRitz test style
// =============================================================================

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

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
// Helpers
// =============================================================================

/**
 * @brief Return a valid, fully-populated OffTarget for reuse across tests.
 *
 * Default is a SpCas9 hit at chr1:100000(+), 1 mismatch, no bulges.
 * Callers may override individual counts to test derived properties.
 */
static OffTarget make_hit(int mm = 1, int bdna = 0, int brna = 0) {
  return OffTarget{
      "chr1",
      100000,
      Strand::Forward,
      "ACGTACGTACGTACGTACGTNGG", // 20 nt guide + NGG placeholder
      "ACGTACGTACGTACGTACGTaGG", // lowercase 'a' = 1 mismatch at pos 20
      mm,
      bdna,
      brna};
}

/**
 * @brief Split a string on a delimiter and return the parts.
 */
static std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> parts;
  std::string part;
  for (char c : s) {
    if (c == delim) {
      parts.push_back(part);
      part.clear();
    } else {
      part += c;
    }
  }
  parts.push_back(part);
  return parts;
}

// =============================================================================
// Strand helpers
// =============================================================================

static void test_to_char() {
  record("to_char(Forward) == '+'", to_char(Strand::Forward) == '+');
  record("to_char(Reverse) == '-'", to_char(Strand::Reverse) == '-');
}

static void test_strand_from_char_valid() {
  record("strand_from_char('+') == Forward",
         strand_from_char('+') == Strand::Forward);
  record("strand_from_char('-') == Reverse",
         strand_from_char('-') == Strand::Reverse);
}

static void test_strand_from_char_invalid() {
  bool threw_on_space = false;
  try {
    strand_from_char(' ');
  } catch (const std::invalid_argument &) {
    threw_on_space = true;
  }
  record("strand_from_char(' ') throws invalid_argument", threw_on_space);

  bool threw_on_dot = false;
  try {
    strand_from_char('.');
  } catch (const std::invalid_argument &) {
    threw_on_dot = true;
  }
  record("strand_from_char('.') throws invalid_argument", threw_on_dot);

  bool threw_on_null = false;
  try {
    strand_from_char('\0');
  } catch (const std::invalid_argument &) {
    threw_on_null = true;
  }
  record("strand_from_char('\\0') throws invalid_argument", threw_on_null);
}

static void test_strand_round_trip() {
  // to_char ∘ strand_from_char is identity on valid characters
  record("round-trip '+'", to_char(strand_from_char('+')) == '+');
  record("round-trip '-'", to_char(strand_from_char('-')) == '-');
}

// =============================================================================
// Construction — valid paths
// =============================================================================

static void test_construction_basic() {
  bool threw = false;
  try {
    auto ot = make_hit();
    (void)ot;
  } catch (...) {
    threw = true;
  }
  record("basic construction does not throw", !threw);
}

static void test_construction_accessors() {
  auto ot = make_hit(2, 1, 0);

  record("chrom accessor", ot.chrom() == "chr1");
  record("pos accessor", ot.pos() == 100000);
  record("strand accessor", ot.strand() == Strand::Forward);
  record("grna accessor", ot.grna() == "ACGTACGTACGTACGTACGTNGG");
  record("target accessor", ot.target() == "ACGTACGTACGTACGTACGTaGG");
  record("mismatches accessor", ot.mismatches() == 2);
  record("bulge_dna accessor", ot.bulge_dna() == 1);
  record("bulge_rna accessor", ot.bulge_rna() == 0);
}

static void test_construction_reverse_strand() {
  OffTarget ot{"chrM", 500, Strand::Reverse, "ACGTNGG", "acgtGGG", 0, 0, 0};
  record("reverse strand stored", ot.strand() == Strand::Reverse);
  record("chrM stored", ot.chrom() == "chrM");
}

static void test_construction_minimum_pos() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", 0, 0, 0};
  } catch (...) {
    threw = true;
  }
  record("pos == 1 is valid", !threw);
}

static void test_construction_zero_edit_distance() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "ACGTNGG", "ACGTGGG", 0, 0, 0};
    (void)ot;
  } catch (...) {
    threw = true;
  }
  record("all-zero edit counts accepted", !threw);
}

static void test_construction_all_bulge_types() {
  // DNA bulge only
  bool ok1 = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A-G", "ACG", 0, 1, 0};
    ok1 = true;
  } catch (...) {
  }
  record("bdna=1 brna=0 accepted", ok1);

  // RNA bulge only
  bool ok2 = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "ACG", "A-G", 0, 0, 1};
    ok2 = true;
  } catch (...) {
  }
  record("bdna=0 brna=1 accepted", ok2);

  // Mixed
  bool ok3 = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A-CG", "AC-G", 0, 1, 1};
    ok3 = true;
  } catch (...) {
  }
  record("bdna=1 brna=1 (mixed) accepted", ok3);
}

// =============================================================================
// Construction — invalid arguments (each guard exercised in isolation)
// =============================================================================

static void test_construction_invalid_empty_chrom() {
  bool threw = false;
  try {
    OffTarget ot{"", 1, Strand::Forward, "A", "A", 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("empty chrom throws invalid_argument", threw);
}

static void test_construction_invalid_pos_zero() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 0, Strand::Forward, "A", "A", 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("pos == 0 throws invalid_argument", threw);
}

static void test_construction_invalid_pos_negative() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", -1, Strand::Forward, "A", "A", 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("pos < 0 throws invalid_argument", threw);
}

static void test_construction_invalid_empty_grna() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "", "A", 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("empty grna throws invalid_argument", threw);
}

static void test_construction_invalid_empty_target() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "", 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("empty target throws invalid_argument", threw);
}

static void test_construction_invalid_negative_mismatches() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", -1, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("mismatches < 0 throws invalid_argument", threw);
}

static void test_construction_invalid_negative_bulge_dna() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", 0, -1, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("bulge_dna < 0 throws invalid_argument", threw);
}

static void test_construction_invalid_negative_bulge_rna() {
  bool threw = false;
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", 0, 0, -1};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("bulge_rna < 0 throws invalid_argument", threw);
}

// =============================================================================
// Derived properties
// =============================================================================

static void test_total_edit_distance() {
  record("0+0+0 == 0", make_hit(0, 0, 0).total_edit_distance() == 0);
  record("1+0+0 == 1", make_hit(1, 0, 0).total_edit_distance() == 1);
  record("0+1+0 == 1", make_hit(0, 1, 0).total_edit_distance() == 1);
  record("0+0+1 == 1", make_hit(0, 0, 1).total_edit_distance() == 1);
  record("2+1+1 == 4", make_hit(2, 1, 1).total_edit_distance() == 4);
  record("3+2+0 == 5", make_hit(3, 2, 0).total_edit_distance() == 5);
}

static void test_has_bulge() {
  record("mm only => !has_bulge", !make_hit(1, 0, 0).has_bulge());
  record("bdna > 0 => has_bulge", make_hit(0, 1, 0).has_bulge());
  record("brna > 0 => has_bulge", make_hit(0, 0, 1).has_bulge());
  record("bdna+brna > 0 => has_bulge", make_hit(0, 1, 1).has_bulge());
  record("zero all => !has_bulge", !make_hit(0, 0, 0).has_bulge());
}

static void test_bulge_type_none() {
  // No bulge → "X"
  record("bulge_type no-bulge == \"X\"", make_hit(1, 0, 0).bulge_type() == "X");
}

static void test_bulge_type_dna() {
  record("bulge_type dna-only == \"DNA\"",
         make_hit(0, 1, 0).bulge_type() == "DNA");
}

static void test_bulge_type_rna() {
  record("bulge_type rna-only == \"RNA\"",
         make_hit(0, 0, 1).bulge_type() == "RNA");
}

static void test_bulge_type_mixed() {
  record("bulge_type mixed == \"DNA,RNA\"",
         make_hit(0, 1, 1).bulge_type() == "DNA,RNA");

  // Larger counts — still mixed
  record("bulge_type (2,3) mixed == \"DNA,RNA\"",
         make_hit(0, 2, 3).bulge_type() == "DNA,RNA");
}

// =============================================================================
// Equality operator
// =============================================================================

static void test_equality_reflexive() {
  auto a = make_hit();
  record("a == a (reflexive)", a == a);
  record("!(a != a)", !(a != a));
}

static void test_equality_same_locus() {
  // Same chrom/pos/strand/target — must be equal regardless of mm/bulge
  auto a = make_hit(1, 0, 0);
  auto b = make_hit(2, 0, 0); // different mismatches, same target string
  record("same locus, different mm => equal", a == b);
  record("same locus, different mm => !!=", !(a != b));
}

static void test_equality_grna_excluded() {
  // Two OffTargets with same locus but different gRNA should still be equal
  // (grna is excluded from equality; target determines the locus)
  OffTarget a{"chr1",
              100000,
              Strand::Forward,
              "ACGTACGTACGTACGTACGTNGG",
              "ACGTACGTACGTACGTACGTaGG",
              1,
              0,
              0};
  OffTarget b{"chr1",
              100000,
              Strand::Forward,
              "GCGTACGTACGTACGTACGTNGG", // different grna
              "ACGTACGTACGTACGTACGTaGG",
              1,
              0,
              0};
  record("different grna, same target => equal", a == b);
}

static void test_equality_different_chrom() {
  OffTarget a{"chr1", 100000, Strand::Forward, "A", "a", 1, 0, 0};
  OffTarget b{"chr2", 100000, Strand::Forward, "A", "a", 1, 0, 0};
  record("different chrom => not equal", a != b);
  record("different chrom => !(a == b)", !(a == b));
}

static void test_equality_different_pos() {
  OffTarget a{"chr1", 100000, Strand::Forward, "A", "a", 0, 0, 0};
  OffTarget b{"chr1", 200000, Strand::Forward, "A", "a", 0, 0, 0};
  record("different pos => not equal", a != b);
}

static void test_equality_different_strand() {
  OffTarget a{"chr1", 100000, Strand::Forward, "A", "a", 0, 0, 0};
  OffTarget b{"chr1", 100000, Strand::Reverse, "A", "a", 0, 0, 0};
  record("different strand => not equal", a != b);
}

static void test_equality_different_target() {
  OffTarget a{"chr1", 100000, Strand::Forward, "A", "AAAA", 0, 0, 0};
  OffTarget b{"chr1", 100000, Strand::Forward, "A", "CCCC", 0, 0, 0};
  record("different target => not equal", a != b);
}

// =============================================================================
// Ordering operator
// =============================================================================

static void test_ordering_irreflexive() {
  auto a = make_hit();
  record("operator< is irreflexive: !(a < a)", !(a < a));
}

static void test_ordering_by_chrom() {
  OffTarget a{"chr1", 1, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget b{"chr2", 1, Strand::Forward, "A", "A", 0, 0, 0};
  record("chr1 < chr2", a < b);
  record("chr2 not < chr1", !(b < a));
}

static void test_ordering_by_pos() {
  OffTarget a{"chr1", 100, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget b{"chr1", 1000, Strand::Forward, "A", "A", 0, 0, 0};
  record("pos 100 < pos 1000", a < b);
  record("pos 1000 not < 100", !(b < a));
}

static void test_ordering_by_strand() {
  // '+' (ASCII 43) < '-' (ASCII 45) — Forward sorts before Reverse
  OffTarget fwd{"chr1", 100, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget rev{"chr1", 100, Strand::Reverse, "A", "A", 0, 0, 0};
  record("Forward < Reverse (same chrom, pos)", fwd < rev);
  record("Reverse not < Forward", !(rev < fwd));
}

static void test_ordering_by_edit_distance() {
  // Same chrom, pos, strand — lower total_edit_distance sorts first
  OffTarget a{"chr1", 100, Strand::Forward, "A", "a", 1, 0, 0}; // TED = 1
  OffTarget b{"chr1", 100, Strand::Forward, "A", "A", 3, 0, 0}; // TED = 3
  // Note: a == b by locus (same target? No — "a" vs "A" differ here)
  // Let them have the same target to test TED tiebreaker
  OffTarget c{"chr1", 100, Strand::Forward, "A", "X", 1, 0, 0}; // TED = 1
  OffTarget d{"chr1", 100, Strand::Forward, "A", "X", 3, 0, 0}; // TED = 3
  record("lower TED sorts first (same locus)", c < d);
  record("higher TED not < lower (same locus)", !(d < c));
}

static void test_ordering_transitivity() {
  OffTarget a{"chr1", 1, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget b{"chr1", 100, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget c{"chr2", 1, Strand::Forward, "A", "A", 0, 0, 0};
  record("transitivity: a < b && b < c => a < c",
         (a < b) && (b < c) && (a < c));
}

static void test_ordering_sort_vector() {
  std::vector<OffTarget> hits;
  hits.push_back({"chr2", 500, Strand::Forward, "A", "A", 0, 0, 0});
  hits.push_back({"chr1", 200, Strand::Reverse, "A", "A", 0, 0, 0});
  hits.push_back({"chr1", 100, Strand::Forward, "A", "A", 0, 0, 0});
  hits.push_back({"chr1", 100, Strand::Forward, "A", "A", 2, 0, 0});

  std::sort(hits.begin(), hits.end());

  record("sorted[0] chrom == chr1", hits[0].chrom() == "chr1");
  record("sorted[0] pos == 100", hits[0].pos() == 100);
  record("sorted[0] strand == fwd", hits[0].strand() == Strand::Forward);
  // hits[0] and hits[1] both land at chr1:100(+); the one with lower TED
  // must sort first among them.
  record("sorted[1] pos == 100 (TED tiebreak within same pos)",
         hits[1].pos() == 100);
  record("sorted[0] TED <= sorted[1] TED (within chr1:100 group)",
         hits[0].total_edit_distance() <= hits[1].total_edit_distance());
  // hits[2] is the chr1:200(-) entry; it comes after both chr1:100 hits.
  record("sorted[2] pos == 200", hits[2].pos() == 200);
  record("sorted[3] chrom == chr2", hits[3].chrom() == "chr2");
}

// =============================================================================
// std::hash
// =============================================================================

static void test_hash_equal_objects_same_hash() {
  auto a = make_hit(1, 0, 0);
  auto b = make_hit(2, 0, 0); // same locus, different mm
  std::hash<OffTarget> h;
  record("equal objects produce same hash", h(a) == h(b));
}

static void test_hash_different_objects_probably_different_hash() {
  OffTarget a{"chr1", 1, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget b{"chr2", 1, Strand::Forward, "A", "A", 0, 0, 0};
  OffTarget c{"chr1", 99, Strand::Forward, "A", "A", 0, 0, 0};
  std::hash<OffTarget> h;
  // Not guaranteed to differ (pigeonhole), but virtually certain for these
  // values
  record("hash(chr1:1) != hash(chr2:1) — expected but not guaranteed",
         h(a) != h(b));
  record("hash(pos:1) != hash(pos:99) — expected but not guaranteed",
         h(a) != h(c));
}

static void test_hash_unordered_set_deduplication() {
  std::unordered_set<OffTarget> seen;

  auto a = make_hit(1);
  auto b = make_hit(3); // same locus as a (same target string)
  OffTarget c{"chr2", 1, Strand::Forward, "A", "A", 0, 0, 0};

  seen.insert(a);
  seen.insert(b); // duplicate — same locus
  seen.insert(c); // different locus

  record("unordered_set deduplicates same-locus hits (size == 2)",
         seen.size() == 2u);
}

static void test_hash_unordered_set_contains_inserted() {
  std::unordered_set<OffTarget> seen;
  auto ot = make_hit();
  seen.insert(ot);
  record("unordered_set::count finds inserted element", seen.count(ot) == 1u);
}

// =============================================================================
// TSV serialization
// =============================================================================

static void test_tsv_header_count() {
  auto hdr = OffTarget::tsv_header();
  record("tsv_header() returns 9 columns", hdr.size() == 9u);
}

static void test_tsv_header_names() {
  auto hdr = OffTarget::tsv_header();
  // Only check that the expected canonical names appear in order.
  // This is the column set documented in to_tsv_row().
  const std::vector<std::string> expected = {
      "chrom",      "pos",       "strand",    "grna",      "target",
      "mismatches", "bulge_dna", "bulge_rna", "bulge_type"};
  bool match = (hdr == expected);
  record("tsv_header() column names match expected order", match);
}

static void test_tsv_row_field_count() {
  auto row = make_hit().to_tsv_row();
  auto fields = split(row, '\t');
  record("to_tsv_row() produces 9 tab-separated fields", fields.size() == 9u);
}

static void test_tsv_row_field_values() {
  OffTarget ot{"chr7", 42, Strand::Reverse, "GCGGNGG", "gCGGaGG", 2, 0, 0};
  auto fields = split(ot.to_tsv_row(), '\t');

  record("TSV field 0 == chrom", fields.size() > 0 && fields[0] == "chr7");
  record("TSV field 1 == pos", fields.size() > 1 && fields[1] == "42");
  record("TSV field 2 == strand", fields.size() > 2 && fields[2] == "-");
  record("TSV field 3 == grna", fields.size() > 3 && fields[3] == "GCGGNGG");
  record("TSV field 4 == target", fields.size() > 4 && fields[4] == "gCGGaGG");
  record("TSV field 5 == mismatches", fields.size() > 5 && fields[5] == "2");
  record("TSV field 6 == bulge_dna", fields.size() > 6 && fields[6] == "0");
  record("TSV field 7 == bulge_rna", fields.size() > 7 && fields[7] == "0");
  record("TSV field 8 == bulge_type", fields.size() > 8 && fields[8] == "X");
}

static void test_tsv_row_bulge_type_strings_in_output() {
  auto check = [](int bdna, int brna,
                  const std::string &expected_type) -> bool {
    OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", 0, bdna, brna};
    auto fields = split(ot.to_tsv_row(), '\t');
    return fields.size() == 9u && fields[8] == expected_type;
  };

  record("TSV bulge_type field 'X'", check(0, 0, "X"));
  record("TSV bulge_type field 'DNA'", check(1, 0, "DNA"));
  record("TSV bulge_type field 'RNA'", check(0, 1, "RNA"));
  record("TSV bulge_type field 'DNA,RNA'", check(1, 1, "DNA,RNA"));
}

static void test_tsv_row_custom_separator() {
  auto row = make_hit().to_tsv_row(',');
  auto fields = split(row, ',');
  record("to_tsv_row(',') produces 9 comma-separated fields",
         fields.size() == 9u);
  record("first CSV field is chrom", !fields.empty() && fields[0] == "chr1");
}

static void test_tsv_row_forward_strand_symbol() {
  OffTarget ot{"chrX", 1, Strand::Forward, "A", "A", 0, 0, 0};
  auto fields = split(ot.to_tsv_row(), '\t');
  record("forward strand serialises as '+'",
         fields.size() > 2 && fields[2] == "+");
}

static void test_tsv_row_reverse_strand_symbol() {
  OffTarget ot{"chrX", 1, Strand::Reverse, "A", "A", 0, 0, 0};
  auto fields = split(ot.to_tsv_row(), '\t');
  record("reverse strand serialises as '-'",
         fields.size() > 2 && fields[2] == "-");
}

static void test_tsv_row_no_trailing_newline() {
  auto row = make_hit().to_tsv_row();
  record("to_tsv_row() produces no trailing newline",
         !row.empty() && row.back() != '\n');
}

static void test_tsv_row_header_row_field_count_matches() {
  auto header_cols = OffTarget::tsv_header();
  auto data_fields = split(make_hit().to_tsv_row(), '\t');
  record("header column count equals data field count",
         header_cols.size() == data_fields.size());
}

// =============================================================================
// locus_string
// =============================================================================

static void test_locus_string_forward() {
  OffTarget ot{"chr1", 123456, Strand::Forward, "A", "A", 0, 0, 0};
  record("locus_string forward == \"chr1:123456(+)\"",
         ot.locus_string() == "chr1:123456(+)");
}

static void test_locus_string_reverse() {
  OffTarget ot{"chrM", 5000, Strand::Reverse, "A", "A", 0, 0, 0};
  record("locus_string reverse == \"chrM:5000(-)\"",
         ot.locus_string() == "chrM:5000(-)");
}

static void test_locus_string_pos_one() {
  OffTarget ot{"chr1", 1, Strand::Forward, "A", "A", 0, 0, 0};
  record("locus_string pos==1 == \"chr1:1(+)\"",
         ot.locus_string() == "chr1:1(+)");
}

// =============================================================================
// Edge cases
// =============================================================================

static void test_edge_max_int32_position() {
  int32_t max_pos = std::numeric_limits<int32_t>::max();
  try {
    OffTarget ot{"chr1", max_pos, Strand::Forward, "A", "A", 0, 0, 0};
    record("INT32_MAX position accepted and stored", ot.pos() == max_pos);
  } catch (...) {
    record("INT32_MAX position accepted and stored", false, "unexpected throw");
  }
}

static void test_edge_long_chromosome_name() {
  std::string long_chrom(256, 'c');
  try {
    OffTarget ot{long_chrom, 1, Strand::Forward, "A", "A", 0, 0, 0};
    record("256-char chrom name accepted", ot.chrom().size() == 256u);
  } catch (...) {
    record("256-char chrom name accepted", false, "unexpected throw");
  }
}

static void test_edge_long_sequences() {
  // 30-nt guide + 3-nt PAM (common for Cas12a)
  std::string grna(33, 'N');
  std::string target(33, 'a');
  try {
    OffTarget ot{"chr1", 1, Strand::Forward, grna, target, 33, 0, 0};
    record("33-nt sequences accepted", ot.grna().size() == 33u);
  } catch (...) {
    record("33-nt sequences accepted", false, "unexpected throw");
  }
}

static void test_edge_exact_match_is_valid() {
  // All zeros for edit distance — a perfect on-target site
  try {
    OffTarget ot{"chr1",
                 1,
                 Strand::Forward,
                 "ACGTACGTACGTACGTACGTNGG",
                 "ACGTACGTACGTACGTACGTGGG",
                 0,
                 0,
                 0};
    record("perfect match (0,0,0) accepted, TED == 0",
           ot.total_edit_distance() == 0);
    record("perfect match has_bulge() == false", !ot.has_bulge());
    record("perfect match bulge_type() == \"X\"", ot.bulge_type() == "X");
  } catch (...) {
    record("perfect match (0,0,0) accepted, TED == 0", false,
           "unexpected throw");
  }
}

static void test_edge_copy_and_compare() {
  // OffTarget is a value type; copies must compare equal to the original
  auto orig = make_hit();
  auto copy = orig; // copy constructor
  record("copy == original", copy == orig);
  record("!(copy != original)", !(copy != orig));

  auto moved = std::move(copy); // move constructor
  record("moved == original", moved == orig);
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_offtarget ===\n\n";

  std::cout << "-- Strand helpers --\n";
  test_to_char();
  test_strand_from_char_valid();
  test_strand_from_char_invalid();
  test_strand_round_trip();

  std::cout << "\n-- Construction (valid) --\n";
  test_construction_basic();
  test_construction_accessors();
  test_construction_reverse_strand();
  test_construction_minimum_pos();
  test_construction_zero_edit_distance();
  test_construction_all_bulge_types();

  std::cout << "\n-- Construction (invalid arguments) --\n";
  test_construction_invalid_empty_chrom();
  test_construction_invalid_pos_zero();
  test_construction_invalid_pos_negative();
  test_construction_invalid_empty_grna();
  test_construction_invalid_empty_target();
  test_construction_invalid_negative_mismatches();
  test_construction_invalid_negative_bulge_dna();
  test_construction_invalid_negative_bulge_rna();

  std::cout << "\n-- Derived properties --\n";
  test_total_edit_distance();
  test_has_bulge();
  test_bulge_type_none();
  test_bulge_type_dna();
  test_bulge_type_rna();
  test_bulge_type_mixed();

  std::cout << "\n-- Equality --\n";
  test_equality_reflexive();
  test_equality_same_locus();
  test_equality_grna_excluded();
  test_equality_different_chrom();
  test_equality_different_pos();
  test_equality_different_strand();
  test_equality_different_target();

  std::cout << "\n-- Ordering --\n";
  test_ordering_irreflexive();
  test_ordering_by_chrom();
  test_ordering_by_pos();
  test_ordering_by_strand();
  test_ordering_by_edit_distance();
  test_ordering_transitivity();
  test_ordering_sort_vector();

  std::cout << "\n-- std::hash --\n";
  test_hash_equal_objects_same_hash();
  test_hash_different_objects_probably_different_hash();
  test_hash_unordered_set_deduplication();
  test_hash_unordered_set_contains_inserted();

  std::cout << "\n-- TSV serialization --\n";
  test_tsv_header_count();
  test_tsv_header_names();
  test_tsv_row_field_count();
  test_tsv_row_field_values();
  test_tsv_row_bulge_type_strings_in_output();
  test_tsv_row_custom_separator();
  test_tsv_row_forward_strand_symbol();
  test_tsv_row_reverse_strand_symbol();
  test_tsv_row_no_trailing_newline();
  test_tsv_row_header_row_field_count_matches();

  std::cout << "\n-- locus_string --\n";
  test_locus_string_forward();
  test_locus_string_reverse();
  test_locus_string_pos_one();

  std::cout << "\n-- Edge cases --\n";
  test_edge_max_int32_position();
  test_edge_long_chromosome_name();
  test_edge_long_sequences();
  test_edge_exact_match_is_valid();
  test_edge_copy_and_compare();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}