/**
 * @file test_nucleotide_encoding.cpp
 * @brief Unit tests for nucleotide_encoding.cpp
 *
 * Covers:
 *  - NucleotideEncoder::encode_genome   – all IUPAC codes + unknown
 *  - NucleotideEncoder::decode_genome   – round-trip against encode_genome
 *  - NucleotideEncoder::encode_pam      – wildcard 'N' vs encode_genome 'N'
 *  - NucleotideEncoder::decode_pam      – round-trip against encode_pam
 *  - NucleotideEncoder::complement      – all IUPAC complements
 *  - pam::reverse_complement            – empty, single, palindrome, mixed
 */

#include "nucleotide_encoding.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

using pam::NucleotideEncoder;
using pam::reverse_complement;

// -----------------------------------------------------------------------------
// Minimal test harness
// -----------------------------------------------------------------------------

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

/**
 * @brief Records one test result; prints a short summary line.
 * @param name   Human-readable test name.
 * @param ok     Whether the test passed.
 * @param detail Optional extra information printed on failure.
 */
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

/** @brief True when encode_genome(c) produces the expected 4-bit pattern. */
static bool genome_encodes_to(char c, uint8_t expected) {
  return NucleotideEncoder::encode_genome(c) == expected;
}

/** @brief True when encode_pam(c) produces the expected 4-bit pattern. */
static bool pam_encodes_to(char c, uint8_t expected) {
  return NucleotideEncoder::encode_pam(c) == expected;
}

// -----------------------------------------------------------------------------
// Test groups
// -----------------------------------------------------------------------------

/** @brief Verify that all canonical IUPAC base characters encode correctly. */
static void test_encode_genome_canonical() {
  record("encode_genome A == 0b0001", genome_encodes_to('A', 0b0001));
  record("encode_genome C == 0b0010", genome_encodes_to('C', 0b0010));
  record("encode_genome G == 0b0100", genome_encodes_to('G', 0b0100));
  record("encode_genome T == 0b1000", genome_encodes_to('T', 0b1000));
  record("encode_genome N == 0b0000", genome_encodes_to('N', 0b0000));
}

/**
 * @brief Verify that all IUPAC ambiguity codes encode to the union of their
 *        constituent bases.
 */
static void test_encode_genome_iupac_ambiguity() {
  // R = A | G
  record("encode_genome R == A|G", genome_encodes_to('R', 0b0101));
  // Y = C | T
  record("encode_genome Y == C|T", genome_encodes_to('Y', 0b1010));
  // S = G | C
  record("encode_genome S == G|C", genome_encodes_to('S', 0b0110));
  // W = A | T
  record("encode_genome W == A|T", genome_encodes_to('W', 0b1001));
  // K = G | T
  record("encode_genome K == G|T", genome_encodes_to('K', 0b1100));
  // M = A | C
  record("encode_genome M == A|C", genome_encodes_to('M', 0b0011));
  // B = C | G | T
  record("encode_genome B == C|G|T", genome_encodes_to('B', 0b1110));
  // D = A | G | T
  record("encode_genome D == A|G|T", genome_encodes_to('D', 0b1101));
  // H = A | C | T
  record("encode_genome H == A|C|T", genome_encodes_to('H', 0b1011));
  // V = A | C | G
  record("encode_genome V == A|C|G", genome_encodes_to('V', 0b0111));
}

/** @brief Unknown characters must not cause UB; they should return 0b0000. */
static void test_encode_genome_unknown() {
  record("encode_genome '?' == 0b0000", genome_encodes_to('?', 0b0000));
  record("encode_genome 'X' == 0b0000", genome_encodes_to('X', 0b0000));
  record("encode_genome ' ' == 0b0000", genome_encodes_to(' ', 0b0000));
  record("encode_genome '0' == 0b0000", genome_encodes_to('0', 0b0000));
}

/**
 * @brief Round-trip: encode_genome followed by decode_genome must return the
 *        original character for all canonical IUPAC codes.
 */
static void test_decode_genome_roundtrip() {
  const std::string codes = "ACGTRYSWKMBDHV";
  for (char c : codes) {
    uint8_t enc = NucleotideEncoder::encode_genome(c);
    char back = NucleotideEncoder::decode_genome(enc);
    record(std::string("genome round-trip ") + c, back == c,
           std::string("got '") + back + "'");
  }
}

/**
 * @brief The PAM encoder must treat 'N' as a wildcard (0b1111) rather than
 *        0b0000 as in the genome encoder.
 */
static void test_encode_pam_n_is_wildcard() {
  record("encode_pam N == 0b1111 (wildcard)", pam_encodes_to('N', 0b1111));
  // All other canonical codes should match their genome encoding
  record("encode_pam A == 0b0001", pam_encodes_to('A', 0b0001));
  record("encode_pam C == 0b0010", pam_encodes_to('C', 0b0010));
  record("encode_pam G == 0b0100", pam_encodes_to('G', 0b0100));
  record("encode_pam T == 0b1000", pam_encodes_to('T', 0b1000));
}

/** @brief Round-trip: encode_pam / decode_pam for ACGT + N. */
static void test_decode_pam_roundtrip() {
  const std::string codes = "ACGTN";
  for (char c : codes) {
    uint8_t enc = NucleotideEncoder::encode_pam(c);
    char back = NucleotideEncoder::decode_pam(enc);
    record(std::string("pam round-trip ") + c, back == c,
           std::string("got '") + back + "'");
  }
}

/** @brief Verify Watson–Crick and IUPAC complement table entries. */
static void test_complement_basic() {
  record("complement A == T", NucleotideEncoder::complement('A') == 'T');
  record("complement T == A", NucleotideEncoder::complement('T') == 'A');
  record("complement C == G", NucleotideEncoder::complement('C') == 'G');
  record("complement G == C", NucleotideEncoder::complement('G') == 'C');
}

/** @brief Verify IUPAC ambiguity complement pairs. */
static void test_complement_iupac_pairs() {
  // R <-> Y
  record("complement R == Y", NucleotideEncoder::complement('R') == 'Y');
  record("complement Y == R", NucleotideEncoder::complement('Y') == 'R');
  // M <-> K
  record("complement M == K", NucleotideEncoder::complement('M') == 'K');
  record("complement K == M", NucleotideEncoder::complement('K') == 'M');
  // H <-> D
  record("complement H == D", NucleotideEncoder::complement('H') == 'D');
  record("complement D == H", NucleotideEncoder::complement('D') == 'H');
  // B <-> V
  record("complement B == V", NucleotideEncoder::complement('B') == 'V');
  record("complement V == B", NucleotideEncoder::complement('V') == 'B');
  // S and W are self-complementary
  record("complement S == S", NucleotideEncoder::complement('S') == 'S');
  record("complement W == W", NucleotideEncoder::complement('W') == 'W');
}

/** @brief Unknown characters should be returned unchanged by complement. */
static void test_complement_unknown_passthrough() {
  record("complement '?' passthrough",
         NucleotideEncoder::complement('?') == '?');
  record("complement 'N' passthrough",
         NucleotideEncoder::complement('N') == 'N');
}

/** @brief reverse_complement of an empty string must be empty. */
static void test_rc_empty() {
  record("rc(\"\") == \"\"", reverse_complement("") == "");
}

/** @brief Single character reverse complements. */
static void test_rc_single_char() {
  record("rc(\"A\") == \"T\"", reverse_complement("A") == "T");
  record("rc(\"T\") == \"A\"", reverse_complement("T") == "A");
  record("rc(\"C\") == \"G\"", reverse_complement("C") == "G");
  record("rc(\"G\") == \"C\"", reverse_complement("G") == "C");
}

/** @brief Palindromic sequence: rc must equal the original. */
static void test_rc_palindrome() {
  // ACGT is self-complementary when reversed: rc(ACGT) = rc(T)rc(G)rc(C)rc(A) =
  // ACGT
  record("rc(\"ACGT\") == \"ACGT\"", reverse_complement("ACGT") == "ACGT");
}

/** @brief Multi-character non-palindromic sequences. */
static void test_rc_general() {
  // rc(AAAA) = TTTT
  record("rc(\"AAAA\") == \"TTTT\"", reverse_complement("AAAA") == "TTTT");
  // rc(ATCG) = CGAT
  record("rc(\"ATCG\") == \"CGAT\"", reverse_complement("ATCG") == "CGAT");
  // rc(GATTACA) = TGTAATC
  record("rc(\"GATTACA\") == \"TGTAATC\"",
         reverse_complement("GATTACA") == "TGTAATC");
}

/** @brief Applying rc twice must recover the original string. */
static void test_rc_involution() {
  const std::vector<std::string> seqs = {"ACGT", "GATTACA", "TTTTAAAA",
                                         "RYMKSWHBVDN", "A"};
  for (const auto &s : seqs) {
    bool ok = reverse_complement(reverse_complement(s)) == s;
    record("rc(rc(" + s + ")) == " + s, ok);
  }
}

/**
 * @brief Bitwise AND of an encoded genome base with its encoded PAM complement
 *        must be non-zero (i.e. each IUPAC code matches at least itself).
 */
static void test_genome_pam_self_match() {
  const std::string codes = "ACGTRYMKSWHBVD";
  for (char c : codes) {
    uint8_t g = NucleotideEncoder::encode_genome(c);
    uint8_t p = NucleotideEncoder::encode_pam(c);
    bool ok = (g & p) != 0;
    // 'N' in genome encodes to 0b0000 so AND with anything is 0 – expected
    if (c == 'N')
      ok = !ok; // flip expectation: N in genome should NOT match
    record(std::string("genome/pam self-match ") + c, ok);
  }
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------

int main() {
  std::cout << "=== test_nucleotide_encoding ===\n\n";

  std::cout << "-- encode_genome canonical --\n";
  test_encode_genome_canonical();

  std::cout << "\n-- encode_genome IUPAC ambiguity --\n";
  test_encode_genome_iupac_ambiguity();

  std::cout << "\n-- encode_genome unknown characters --\n";
  test_encode_genome_unknown();

  std::cout << "\n-- decode_genome round-trip --\n";
  test_decode_genome_roundtrip();

  std::cout << "\n-- encode_pam N is wildcard --\n";
  test_encode_pam_n_is_wildcard();

  std::cout << "\n-- decode_pam round-trip --\n";
  test_decode_pam_roundtrip();

  std::cout << "\n-- complement basic --\n";
  test_complement_basic();

  std::cout << "\n-- complement IUPAC pairs --\n";
  test_complement_iupac_pairs();

  std::cout << "\n-- complement unknown passthrough --\n";
  test_complement_unknown_passthrough();

  std::cout << "\n-- reverse_complement edge cases --\n";
  test_rc_empty();
  test_rc_single_char();
  test_rc_palindrome();
  test_rc_general();
  test_rc_involution();

  std::cout << "\n-- genome/pam self-match --\n";
  test_genome_pam_self_match();

  std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}
