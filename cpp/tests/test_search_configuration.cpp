/**
 * @file test_search_configuration.cpp
 * @brief Unit tests for the search configuration layer
 *        (search_configuration.hpp / search_configuration.cpp).
 *
 * Scope note: SearchConfiguration intentionally does NOT carry the PAM or the
 * guide length. The PAM is applied at indexing time and its geometry is stored
 * in the .bin header; the edit-budget-vs-guide-length check is therefore
 * deferred to the index loader. These tests assert that deferral as a positive
 * property (a large edit budget constructs fine standalone) rather than testing
 * a check that no longer lives here.
 *
 * Test surface:
 *   - OutputFormat round-trip (to_string / output_format_from_string)
 *   - OutputMode round-trip   (to_string / output_mode_from_string)
 *   - SearchConfiguration construction (valid, invalid, boundary)
 *   - Accessors
 *   - Derived: max_bulges_total, max_total_edits, bulges_disabled
 *   - Default output_format (Tsv) and output_mode (Both)
 *   - Convenience: write_targets() / write_profile()
 *   - Deferred guide-length check (large budget accepted standalone)
 *   - Copy/move value semantics
 *
 * Uses the same lightweight record()/g_* harness as the other CRISPRitz C++
 * tests so it registers with add_crispritz_test without an external framework.
 */

#include "search_configuration.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using crispritz::output_format_from_string;
using crispritz::output_mode_from_string;
using crispritz::OutputFormat;
using crispritz::OutputMode;
using crispritz::SearchConfiguration;
using crispritz::to_string;

// =============================================================================
// Minimal test harness
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

/** @brief A default valid config: 4 mismatches, no bulges, single thread. */
static SearchConfiguration make_config() {
  return SearchConfiguration{/*mm=*/4, /*bdna=*/0, /*brna=*/0, /*threads=*/1};
}

// =============================================================================
// OutputFormat
// =============================================================================

static void test_output_format_to_string() {
  record("to_string(Tsv) == \"tsv\"", to_string(OutputFormat::Tsv) == "tsv");
  record("to_string(Targets) == \"targets\"",
         to_string(OutputFormat::Targets) == "targets");
}

static void test_output_format_from_string_valid() {
  record("from_string(\"tsv\") == Tsv",
         output_format_from_string("tsv") == OutputFormat::Tsv);
  record("from_string(\"targets\") == Targets",
         output_format_from_string("targets") == OutputFormat::Targets);
}

static void test_output_format_from_string_invalid() {
  bool threw = false;
  try {
    (void)output_format_from_string("json");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("from_string(\"json\") throws invalid_argument", threw);

  bool threw_empty = false;
  try {
    (void)output_format_from_string("");
  } catch (const std::invalid_argument &) {
    threw_empty = true;
  }
  record("from_string(\"\") throws invalid_argument", threw_empty);

  bool threw_case = false;
  try {
    (void)output_format_from_string("TSV");
  } // case-sensitive
  catch (const std::invalid_argument &) {
    threw_case = true;
  }
  record("from_string(\"TSV\") throws (case-sensitive)", threw_case);
}

static void test_output_format_round_trip() {
  record("round-trip Tsv", output_format_from_string(to_string(
                               OutputFormat::Tsv)) == OutputFormat::Tsv);
  record("round-trip Targets",
         output_format_from_string(to_string(OutputFormat::Targets)) ==
             OutputFormat::Targets);
}

// =============================================================================
// OutputMode
// =============================================================================

static void test_output_mode_to_string() {
  record("to_string(TargetsOnly) == \"targets\"",
         to_string(OutputMode::TargetsOnly) == "targets");
  record("to_string(ProfileOnly) == \"profile\"",
         to_string(OutputMode::ProfileOnly) == "profile");
  record("to_string(Both) == \"both\"", to_string(OutputMode::Both) == "both");
}

static void test_output_mode_from_string_valid() {
  record("from_string(\"targets\") == TargetsOnly",
         output_mode_from_string("targets") == OutputMode::TargetsOnly);
  record("from_string(\"profile\") == ProfileOnly",
         output_mode_from_string("profile") == OutputMode::ProfileOnly);
  record("from_string(\"both\") == Both",
         output_mode_from_string("both") == OutputMode::Both);
}

static void test_output_mode_from_string_invalid() {
  bool threw = false;
  try {
    (void)output_mode_from_string("all");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("from_string(\"all\") throws invalid_argument", threw);

  bool threw_empty = false;
  try {
    (void)output_mode_from_string("");
  } catch (const std::invalid_argument &) {
    threw_empty = true;
  }
  record("from_string(\"\") throws invalid_argument", threw_empty);

  bool threw_case = false;
  try {
    (void)output_mode_from_string("Both");
  } catch (const std::invalid_argument &) {
    threw_case = true;
  }
  record("from_string(\"Both\") throws (case-sensitive)", threw_case);
}

static void test_output_mode_round_trip() {
  record("round-trip TargetsOnly",
         output_mode_from_string(to_string(OutputMode::TargetsOnly)) ==
             OutputMode::TargetsOnly);
  record("round-trip ProfileOnly",
         output_mode_from_string(to_string(OutputMode::ProfileOnly)) ==
             OutputMode::ProfileOnly);
  record("round-trip Both", output_mode_from_string(to_string(
                                OutputMode::Both)) == OutputMode::Both);
}

// =============================================================================
// SearchConfiguration — valid construction
// =============================================================================

static void test_config_valid() {
  bool threw = false;
  try {
    auto c = make_config();
    (void)c;
  } catch (...) {
    threw = true;
  }
  record("basic config constructs", !threw);
}

static void test_config_accessors() {
  SearchConfiguration c{3, 1, 1, 8, OutputFormat::Targets};
  record("config max_mismatches == 3", c.max_mismatches() == 3);
  record("config max_bulges_dna == 1", c.max_bulges_dna() == 1);
  record("config max_bulges_rna == 1", c.max_bulges_rna() == 1);
  record("config threads == 8", c.threads() == 8);
  record("config output_format == Targets",
         c.output_format() == OutputFormat::Targets);
}

static void test_config_default_output_format() {
  SearchConfiguration c{4, 0, 0, 1};
  record("default output_format == Tsv",
         c.output_format() == OutputFormat::Tsv);
}

static void test_config_default_output_mode() {
  // Default output_mode must be Both so existing callers keep full output.
  SearchConfiguration c{4, 0, 0, 1};
  record("default output_mode == Both", c.output_mode() == OutputMode::Both);
}

static void test_config_explicit_output_mode() {
  SearchConfiguration targets_only{
      4, 0, 0, 1, OutputFormat::Tsv, OutputMode::TargetsOnly};
  record("explicit TargetsOnly stored",
         targets_only.output_mode() == OutputMode::TargetsOnly);

  SearchConfiguration profile_only{
      4, 0, 0, 1, OutputFormat::Tsv, OutputMode::ProfileOnly};
  record("explicit ProfileOnly stored",
         profile_only.output_mode() == OutputMode::ProfileOnly);

  SearchConfiguration both{4, 0, 0, 1, OutputFormat::Tsv, OutputMode::Both};
  record("explicit Both stored", both.output_mode() == OutputMode::Both);
}

// =============================================================================
// SearchConfiguration — write_targets / write_profile convenience helpers
// =============================================================================

static void test_write_targets_helper() {
  SearchConfiguration t{4, 0, 0, 1, OutputFormat::Tsv, OutputMode::TargetsOnly};
  record("write_targets(): TargetsOnly → true", t.write_targets());
  record("write_profile(): TargetsOnly → false", !t.write_profile());

  SearchConfiguration p{4, 0, 0, 1, OutputFormat::Tsv, OutputMode::ProfileOnly};
  record("write_targets(): ProfileOnly → false", !p.write_targets());
  record("write_profile(): ProfileOnly → true", p.write_profile());

  SearchConfiguration b{4, 0, 0, 1, OutputFormat::Tsv, OutputMode::Both};
  record("write_targets(): Both → true", b.write_targets());
  record("write_profile(): Both → true", b.write_profile());
}

// =============================================================================
// SearchConfiguration — derived properties
// =============================================================================

static void test_config_derived_bulges_total() {
  SearchConfiguration c{2, 2, 3, 1};
  record("max_bulges_total == 5", c.max_bulges_total() == 5);
}

static void test_config_derived_total_edits() {
  SearchConfiguration c{2, 2, 3, 1};
  record("max_total_edits == 7 (2+2+3)", c.max_total_edits() == 7);
}

static void test_config_derived_bulges_disabled() {
  SearchConfiguration no_bulge{4, 0, 0, 1};
  record("no bulges => bulges_disabled() true", no_bulge.bulges_disabled());

  SearchConfiguration with_dna{4, 1, 0, 1};
  record("dna bulge => bulges_disabled() false", !with_dna.bulges_disabled());

  SearchConfiguration with_rna{4, 0, 1, 1};
  record("rna bulge => bulges_disabled() false", !with_rna.bulges_disabled());
}

// =============================================================================
// SearchConfiguration — invalid construction
// =============================================================================

static void test_config_invalid_negative_mismatches() {
  bool threw = false;
  try {
    SearchConfiguration c{-1, 0, 0, 1};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("max_mismatches < 0 throws", threw);
}

static void test_config_invalid_negative_bulges_dna() {
  bool threw = false;
  try {
    SearchConfiguration c{0, -1, 0, 1};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("max_bulges_dna < 0 throws", threw);
}

static void test_config_invalid_negative_bulges_rna() {
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, -1, 1};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("max_bulges_rna < 0 throws", threw);
}

static void test_config_invalid_threads_zero() {
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, 0, 0};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("threads == 0 throws", threw);
}

static void test_config_invalid_threads_negative() {
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, 0, -4};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("threads < 0 throws", threw);
}

// =============================================================================
// SearchConfiguration — boundary values
// =============================================================================

static void test_config_boundary_zero_edits() {
  // Exact-match-only search: all budgets zero, single thread.
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, 0, 1};
    (void)c;
  } catch (...) {
    threw = true;
  }
  record("zero edit budget constructs", !threw);
}

static void test_config_boundary_minimum_threads() {
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, 0, 1};
    (void)c;
  } catch (...) {
    threw = true;
  }
  record("threads == 1 constructs", !threw);
}

static void test_config_boundary_high_threads() {
  bool threw = false;
  try {
    SearchConfiguration c{0, 0, 0, 256};
    (void)c;
  } catch (...) {
    threw = true;
  }
  record("threads == 256 constructs", !threw);
}

// =============================================================================
// Deferred guide-length check
// =============================================================================

static void test_config_large_budget_accepted_standalone() {
  // The edit-budget-vs-guide-length check is the index loader's
  // responsibility.
  bool threw = false;
  try {
    SearchConfiguration c{1000, 500, 500, 1};
    record("large budget: max_total_edits == 2000",
           c.max_total_edits() == 2000);
  } catch (...) {
    threw = true;
  }
  record("large edit budget accepted standalone (guide-length check deferred)",
         !threw);
}

// =============================================================================
// Value semantics
// =============================================================================

static void test_config_copy_semantics() {
  SearchConfiguration orig{
      3, 1, 2, 4, OutputFormat::Targets, OutputMode::ProfileOnly};
  auto copy = orig;
  record("copy preserves max_mismatches",
         copy.max_mismatches() == orig.max_mismatches());
  record("copy preserves max_bulges_total",
         copy.max_bulges_total() == orig.max_bulges_total());
  record("copy preserves output_format",
         copy.output_format() == orig.output_format());
  record("copy preserves output_mode",
         copy.output_mode() == orig.output_mode());
}

static void test_config_move_semantics() {
  SearchConfiguration orig{3, 1, 2, 4, OutputFormat::Tsv, OutputMode::Both};
  const int total = orig.max_total_edits();
  const OutputMode mode = orig.output_mode();
  auto moved = std::move(orig);
  record("moved preserves max_total_edits", moved.max_total_edits() == total);
  record("moved preserves threads", moved.threads() == 4);
  record("moved preserves output_mode", moved.output_mode() == mode);
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_search_configuration ===\n\n";

  std::cout << "-- OutputFormat --\n";
  test_output_format_to_string();
  test_output_format_from_string_valid();
  test_output_format_from_string_invalid();
  test_output_format_round_trip();

  std::cout << "\n-- OutputMode --\n";
  test_output_mode_to_string();
  test_output_mode_from_string_valid();
  test_output_mode_from_string_invalid();
  test_output_mode_round_trip();

  std::cout << "\n-- SearchConfiguration (valid) --\n";
  test_config_valid();
  test_config_accessors();
  test_config_default_output_format();
  test_config_default_output_mode();
  test_config_explicit_output_mode();

  std::cout << "\n-- SearchConfiguration (write_targets / write_profile) --\n";
  test_write_targets_helper();

  std::cout << "\n-- SearchConfiguration (derived) --\n";
  test_config_derived_bulges_total();
  test_config_derived_total_edits();
  test_config_derived_bulges_disabled();

  std::cout << "\n-- SearchConfiguration (invalid) --\n";
  test_config_invalid_negative_mismatches();
  test_config_invalid_negative_bulges_dna();
  test_config_invalid_negative_bulges_rna();
  test_config_invalid_threads_zero();
  test_config_invalid_threads_negative();

  std::cout << "\n-- SearchConfiguration (boundary) --\n";
  test_config_boundary_zero_edits();
  test_config_boundary_minimum_threads();
  test_config_boundary_high_threads();

  std::cout << "\n-- Deferred guide-length check --\n";
  test_config_large_budget_accepted_standalone();

  std::cout << "\n-- Value semantics --\n";
  test_config_copy_semantics();
  test_config_move_semantics();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}