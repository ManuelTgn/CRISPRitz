/**
 * @file test_tst_utils.cpp
 * @brief Unit tests for tst_utils.hpp / tst_utils.cpp
 *
 * Covers:
 *  - pack_nibbles / high_nibble / low_nibble round-trips
 *  - SENTINEL_NIBBLE and NULL_CHILD_NIBBLE constant values
 *  - iupac::encode_genome alias agreement with NucleotideEncoder::encode_genome
 *  - iupac::encode_pam alias agreement with NucleotideEncoder::encode_pam
 *  - iupac::complement alias agreement with NucleotideEncoder::complement
 *  - iupac::matches helper semantics
 *  - crispritz::reverse_complement alias correctness
 */

#include "nucleotide_encoding.hpp"
#include "tst_utils.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace crispritz;

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
// pack_nibbles / high_nibble / low_nibble
// -----------------------------------------------------------------------------

/** @brief pack_nibbles(h, l)[7:4] == h, [3:0] == l for all 4-bit combos. */
static void test_pack_nibbles_all_values() {
  bool ok = true;
  for (uint8_t h = 0; h < 16; ++h) {
    for (uint8_t l = 0; l < 16; ++l) {
      uint8_t packed = pack_nibbles(h, l);
      if (high_nibble(packed) != h || low_nibble(packed) != l) {
        ok = false;
        std::cout << "    h=" << (int)h << " l=" << (int)l
                  << " packed=" << (int)packed
                  << " hi=" << (int)high_nibble(packed)
                  << " lo=" << (int)low_nibble(packed) << "\n";
      }
    }
  }
  record("pack/unpack round-trip for all 256 pairs", ok);
}

/** @brief high_nibble of a byte whose upper half is 0x0F must be 0x0F. */
static void test_high_nibble_extraction() {
  record("high_nibble(0xF0) == 0x0F", high_nibble(0xF0) == 0x0F);
  record("high_nibble(0xA0) == 0x0A", high_nibble(0xA0) == 0x0A);
  record("high_nibble(0x00) == 0x00", high_nibble(0x00) == 0x00);
}

/** @brief low_nibble of a byte whose lower half is 0x0F must be 0x0F. */
static void test_low_nibble_extraction() {
  record("low_nibble(0x0F) == 0x0F", low_nibble(0x0F) == 0x0F);
  record("low_nibble(0x0A) == 0x0A", low_nibble(0x0A) == 0x0A);
  record("low_nibble(0x00) == 0x00", low_nibble(0x00) == 0x00);
}

// -----------------------------------------------------------------------------
// Sentinel / null-child constants
// -----------------------------------------------------------------------------

/** @brief SENTINEL_NIBBLE must be 0b1111 (0x0F). */
static void test_sentinel_nibble_value() {
  record("SENTINEL_NIBBLE == 0x0F",
         SENTINEL_NIBBLE == static_cast<uint8_t>(0x0F));
}

/** @brief NULL_CHILD_NIBBLE must be 0b0000 (0x00). */
static void test_null_child_nibble_value() {
  record("NULL_CHILD_NIBBLE == 0x00",
         NULL_CHILD_NIBBLE == static_cast<uint8_t>(0x00));
}

/** @brief Sentinel and null must be distinct. */
static void test_sentinel_vs_null_distinct() {
  record("SENTINEL_NIBBLE != NULL_CHILD_NIBBLE",
         SENTINEL_NIBBLE != NULL_CHILD_NIBBLE);
}

// -----------------------------------------------------------------------------
// iupac:: aliases vs NucleotideEncoder direct calls
// -----------------------------------------------------------------------------

/**
 * @brief iupac::encode_genome must agree with NucleotideEncoder::encode_genome
 *        for all canonical IUPAC characters.
 */
static void test_iupac_encode_genome_alias() {
  const std::string codes = "ACGTNRYMKSWHBVD";
  bool ok = true;
  for (char c : codes) {
    uint8_t via_alias = iupac::encode_genome(c);
    uint8_t via_direct = pam::NucleotideEncoder::encode_genome(c);
    if (via_alias != via_direct) {
      ok = false;
      std::cout << "    mismatch for '" << c << "': alias=" << (int)via_alias
                << " direct=" << (int)via_direct << "\n";
    }
  }
  record("iupac::encode_genome alias agrees with NucleotideEncoder", ok);
}

/**
 * @brief iupac::encode_pam must agree with NucleotideEncoder::encode_genome
 *        for ACGTN (most common PAM characters).
 */
static void test_iupac_encode_pam_alias() {
  const std::string codes = "ACGTN";
  bool ok = true;
  for (char c : codes) {
    uint8_t via_alias = iupac::encode_pam(c);
    uint8_t via_direct = pam::NucleotideEncoder::encode_genome(c);
    if (via_alias != via_direct) {
      ok = false;
      std::cout << "    mismatch for '" << c << "'\n";
    }
  }
  record("iupac::encode_pam alias agrees with NucleotideEncoder", ok);
}

/**
 * @brief iupac::complement must agree with NucleotideEncoder::complement
 *        for ACGT and IUPAC pairs.
 */
static void test_iupac_complement_alias() {
  const std::string codes = "ACGTRYSWMKHDBV";
  bool ok = true;
  for (char c : codes) {
    char via_alias = iupac::complement(c);
    char via_direct = pam::NucleotideEncoder::complement(c);
    if (via_alias != via_direct) {
      ok = false;
      std::cout << "    mismatch for '" << c << "'\n";
    }
  }
  record("iupac::complement alias agrees with NucleotideEncoder", ok);
}

// -----------------------------------------------------------------------------
// iupac::matches
// -----------------------------------------------------------------------------

/** @brief A genome base should match itself (AND of same code is non-zero). */
static void test_iupac_matches_self() {
  const std::string bases = "ACGT";
  bool ok = true;
  for (char c : bases) {
    uint8_t enc = iupac::encode_genome(c);
    if (!iupac::matches(enc, enc)) {
      ok = false;
      std::cout << "    " << c << " does not match itself\n";
    }
  }
  record("iupac::matches: each canonical base matches itself", ok);
}

/** @brief A base should NOT match its complement (e.g. A vs T). */
static void test_iupac_matches_complement_no_match() {
  // A (0001) & T (1000) == 0
  record("iupac::matches(A, T) == false",
         !iupac::matches(iupac::encode_genome('A'), iupac::encode_genome('T')));
  record("iupac::matches(C, G) == false",
         !iupac::matches(iupac::encode_genome('C'), iupac::encode_genome('G')));
}

/** @brief R (A|G) must match both A and G individually. */
static void test_iupac_matches_ambiguity() {
  uint8_t R = iupac::encode_genome('R'); // 0b0101
  uint8_t A = iupac::encode_genome('A'); // 0b0001
  uint8_t G = iupac::encode_genome('G'); // 0b0100
  uint8_t C = iupac::encode_genome('C'); // 0b0010
  record("iupac::matches(R, A) == true", iupac::matches(R, A));
  record("iupac::matches(R, G) == true", iupac::matches(R, G));
  record("iupac::matches(R, C) == false", !iupac::matches(R, C));
}

/** @brief N in genome (0b0000) matches nothing (AND is always 0). */
static void test_iupac_matches_n_genome_matches_nothing() {
  uint8_t N_genome = iupac::encode_genome('N'); // 0b0000
  const std::string bases = "ACGT";
  bool ok = true;
  for (char c : bases) {
    if (iupac::matches(N_genome, iupac::encode_genome(c))) {
      ok = false;
      std::cout << "    N in genome erroneously matched '" << c << "'\n";
    }
  }
  record("N in genome does not match any canonical base", ok);
}

// -----------------------------------------------------------------------------
// crispritz::reverse_complement alias
// -----------------------------------------------------------------------------

/** @brief The alias must forward to pam::reverse_complement correctly. */
static void test_rc_alias_basic() {
  record("crispritz::rc(\"ACGT\") == \"ACGT\"",
         reverse_complement("ACGT") == "ACGT");
  record("crispritz::rc(\"AAAA\") == \"TTTT\"",
         reverse_complement("AAAA") == "TTTT");
  record("crispritz::rc(\"GATTACA\") == \"TGTAATC\"",
         reverse_complement("GATTACA") == "TGTAATC");
}

/** @brief Applying the alias twice recovers the original sequence. */
static void test_rc_alias_involution() {
  const std::vector<std::string> seqs = {"ACGT", "GATTACA", "A", "RYMKSWHBVDN"};
  for (const auto &s : seqs) {
    record("crispritz::rc(rc(" + s + ")) == " + s,
           reverse_complement(reverse_complement(s)) == s);
  }
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------

int main() {
  std::cout << "=== test_tst_utils ===\n\n";

  std::cout << "-- pack_nibbles / high_nibble / low_nibble --\n";
  test_pack_nibbles_all_values();
  test_high_nibble_extraction();
  test_low_nibble_extraction();

  std::cout << "\n-- sentinel / null constants --\n";
  test_sentinel_nibble_value();
  test_null_child_nibble_value();
  test_sentinel_vs_null_distinct();

  std::cout << "\n-- iupac:: aliases vs NucleotideEncoder --\n";
  test_iupac_encode_genome_alias();
  test_iupac_encode_pam_alias();
  test_iupac_complement_alias();

  std::cout << "\n-- iupac::matches --\n";
  test_iupac_matches_self();
  test_iupac_matches_complement_no_match();
  test_iupac_matches_ambiguity();
  test_iupac_matches_n_genome_matches_nothing();

  std::cout << "\n-- crispritz::reverse_complement alias --\n";
  test_rc_alias_basic();
  test_rc_alias_involution();

  std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}
