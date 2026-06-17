/**
 * @file test_profile_data.cpp
 * @brief Unit tests for ProfileAccumulator and GuideProfile
 *        (profile_data.hpp / profile_data.cpp).
 *
 * ## Alignment conventions used throughout
 *
 * All tests use a 5-bp guide body "ACGTA" with a 3-bp PAM "NGG".
 * Unless otherwise stated, pam_at_start=false, so the alignment string
 * layout is:
 *
 *   [5 body columns][3 PAM columns]   total aln_len = 8 (no bulge)
 *                                                   = 9 (1 DNA bulge)
 *
 * For DNA bulge: aln_len = guide_len + bulge_dna + pam_len.
 * For RNA bulge: aln_len = guide_len            + pam_len  (body stays 5 cols).
 *
 * Gap conventions in grna / target strings:
 *   grna[i]   == '-'  →  DNA bulge (target has extra base, guide has gap)
 *   target[i] == '-'  →  RNA bulge (guide has extra base, target has gap)
 *
 * ## Which counter belongs to which output file
 *
 *  File 1 (.profile.xls)          : pos_mm_count (MM-only mismatches),
 *                                    offt_by_mm, ont_count
 *  File 2 (.extended_profile.xls) : ext_mm_nuc_pos, ext_total_by_mm,
 *                                    ext_dna_by_mm_pos, ext_rna_by_mm_pos
 *  File 3 (.profile_dna.xls)      : pos_bulge_dna, pos_mm_in_dna,
 *                                    offt_dna, ont_count_dna
 *  File 4 (.profile_rna.xls)      : pos_bulge_rna, pos_mm_in_rna,
 *                                    offt_rna, ont_count_rna
 *  File 5 (.profile_complete.xls) : pos_mm_complete, offt_complete_by_mm,
 *                                    ont_count_complete
 *
 * ## Test coverage
 *
 *   - Construction: valid geometry, invalid argument rejection
 *   - Build on empty accumulator: geometry fields, array sizes, all-zero counts
 *   - MM-only push: on-target, single-position mismatch, two-position mismatch
 *   - Nucleotide routing: A/C/G/T each land in the correct ext_mm_nuc_pos slot
 *   - Ambiguous 'N': not counted as a mismatch
 *   - DNA-bulge push: positional bulge count, no pos_mm_count pollution,
 *                     bulge channel bucketing, mismatch within bulge hit
 *   - RNA-bulge push: symmetric coverage to DNA
 *   - PAM handling: pam_at_start=false (body first) and pam_at_start=true (PAM
 * first)
 *   - Complete channel: aggregates all channels; ont_count_complete = any 0-MM
 * hit
 *   - Multi-hit accumulation: counts sum correctly across push() calls
 *   - build() idempotency: two calls return identical profiles
 *   - reset(): clears all counters, preserves geometry; fresh build is all-zero
 *   - GuideProfile value semantics: copy and move
 *   - ext_total consistency: ext_total_by_mm[mm][pos] == Σ
 * ext_mm_nuc_pos[mm][nuc][pos]
 *   - mm clamping: hit with mm > max_mm is silently clamped, no crash
 */

#include "offtarget.hpp"
#include "profile_data.hpp"

#include <iostream>
#include <numeric> // std::accumulate
#include <stdexcept>
#include <string>
#include <vector>

using crispritz::GuideProfile;
using crispritz::OffTarget;
using crispritz::ProfileAccumulator;
using crispritz::Strand;

// =============================================================================
// Minimal test harness (same pattern as the other CRISPRitz test suites)
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
// Fixtures
// =============================================================================

// Standard accumulator used in most tests.
//   guide body = "ACGTA" (5 bp), PAM = "NGG" (3 bp), PAM at end, max_mm=2,
//   max_bulge_dna=1, max_bulge_rna=1
static ProfileAccumulator make_acc(int max_mm = 2, int max_bd = 1,
                                   int max_br = 1, bool pam_start = false) {
  return ProfileAccumulator{"ACGTA", 5, 3, max_mm, max_bd, max_br, pam_start};
}

// ── OffTarget factory helpers ───────────────────────────────────────────────
// Alignment format (pam_at_end): [5 body cols][3 PAM cols]  aln_len = 8

/// Perfect match (on-target, 0 MM, 0 bulge).
static OffTarget ot_exact() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "ACGTANGG", 0, 0, 0};
}

/// 1 MM at body position 0 (guide A → target C).
static OffTarget ot_1mm_pos0_C() {
  // grna[0]='A', target[0]='c' (lowercase mismatch convention).
  // toupper: A vs C → mismatch; nuc_index('C') = 1.
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "cCGTANGG", 1, 0, 0};
}

/// 1 MM at body position 2 (guide G → target A).
static OffTarget ot_1mm_pos2_A() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "ACaTANGG", 1, 0, 0};
}

/// 1 MM at body position 4 (guide A → target T).
static OffTarget ot_1mm_pos4_T() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "ACGTtNGG", 1, 0, 0};
}

/// 2 MM: body pos 1 (C→G) and body pos 3 (T→A).
static OffTarget ot_2mm_pos1G_pos3A() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "AgGaANGG", 2, 0, 0};
}

/// MM at pos 2 with each possible off-target nucleotide.
static OffTarget ot_1mm_pos2_nuc(char nuc) {
  // guide[2]='G'; target[2] = nuc (lowercase for legibility).
  std::string target = "AC_TANGG";
  target[2] = static_cast<char>(std::tolower(static_cast<unsigned char>(nuc)));
  return {"chr1", 100, Strand::Forward, "ACGTANGG", target, 1, 0, 0};
}

/// Hit where target has 'N' at body pos 2 (ambiguous — must not count as MM).
static OffTarget ot_N_at_pos2() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "ACnTANGG", 0, 0, 0};
}

// ── DNA-bulge helpers (aln_len = 9, body_end = 6) ──────────────────────────
// Layout: [6 body cols][3 PAM cols]
// grna has '-' where guide has a gap (extra target base).

/// DNA bulge at body pos 1, 0 MM.
/// Trace: A(→1) –(DNA@1) C(→2) G(→3) T(→4) A(→5)
static OffTarget ot_dna_bulge_pos1_0mm() {
  return {"chr1", 100, Strand::Forward, "A-CGTANGG", "ATCGTANGG", 0, 1, 0};
}

/// DNA bulge at body pos 1, 1 MM at body pos 3 (T→A).
/// Trace: A(→1) –(DNA@1,body_pos stays 1) C(→2) G(→3) T/a mismatch(→4) A(→5)
static OffTarget ot_dna_bulge_pos1_1mm_pos3() {
  return {"chr1", 100, Strand::Forward, "A-CGTANGG", "ATCGaANGG", 1, 1, 0};
}

// ── RNA-bulge helpers (aln_len = 8, body_end = 5) ──────────────────────────
// Layout: [5 body cols][3 PAM cols]
// target has '-' where guide has an extra base.

/// RNA bulge at body pos 2, 0 MM.
/// Trace: A(→1) C(→2) G/-(RNA@2)(→3) T(→4) A(→5)
static OffTarget ot_rna_bulge_pos2_0mm() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "AC-TANGG", 0, 0, 1};
}

/// RNA bulge at body pos 1, 1 MM at body pos 3 (T→A).
/// Trace: A(→1) C/-(RNA@1)(→2) G(→3) T/a mismatch(→4) A(→5)
static OffTarget ot_rna_bulge_pos1_1mm_pos3() {
  return {"chr1", 100, Strand::Forward, "ACGTANGG", "A-GaANGG", 1, 0, 1};
}

// ── PAM-at-start helpers (aln_len = 8, body_start = 3, body_end = 8) ────────

/// PAM-at-start, perfect match.
static OffTarget ot_pam_start_exact() {
  return {"chr1", 100, Strand::Forward, "NGGACGTA", "NGGACGTA", 0, 0, 0};
}

/// PAM-at-start, 1 MM at body pos 2 (guide G → target A).
/// body_start=3; i=5 (=3+2) is the mismatch column.
static OffTarget ot_pam_start_1mm_pos2() {
  return {"chr1", 100, Strand::Forward, "NGGACGTA", "NGGACaTA", 1, 0, 0};
}

// =============================================================================
// 1. Construction
// =============================================================================

static void test_construction_valid() {
  bool threw = false;
  try {
    auto a = make_acc();
    (void)a;
  } catch (...) {
    threw = true;
  }
  record("construction: valid params succeed", !threw);
}

static void test_construction_invalid_guide_len_zero() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", 0, 3, 2, 1, 1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: guide_len=0 throws", threw);
}

static void test_construction_invalid_guide_len_negative() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", -1, 3, 2, 1, 1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: guide_len<0 throws", threw);
}

static void test_construction_invalid_pam_len_negative() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", 5, -1, 2, 1, 1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: pam_len<0 throws", threw);
}

static void test_construction_invalid_max_mm_negative() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", 5, 3, -1, 1, 1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: max_mm<0 throws", threw);
}

static void test_construction_invalid_max_bulge_dna_negative() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", 5, 3, 2, -1, 1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: max_bulge_dna<0 throws", threw);
}

static void test_construction_invalid_max_bulge_rna_negative() {
  bool threw = false;
  try {
    ProfileAccumulator{"ACGTA", 5, 3, 2, 1, -1, false};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("construction: max_bulge_rna<0 throws", threw);
}

static void test_construction_boundary_zero_mm_no_bulge() {
  // Exact-match-only search: max_mm=0, no bulges allowed.
  bool threw = false;
  try {
    ProfileAccumulator{"A", 1, 0, 0, 0, 0, false};
  } catch (...) {
    threw = true;
  }
  record("construction: min geometry (guide_len=1,pam_len=0,max_mm=0) succeeds",
         !threw);
}

// =============================================================================
// 2. build() on empty accumulator — geometry and sizes
// =============================================================================

static void test_build_empty_geometry() {
  auto p = make_acc(/*max_mm=*/3, /*max_bd=*/2, /*max_br=*/2).build();

  record("build/empty: guide preserved", p.guide == "ACGTA");
  record("build/empty: guide_len == 5", p.guide_len == 5);
  record("build/empty: pam_len == 3", p.pam_len == 3);
  record("build/empty: pam_at_start == false", !p.pam_at_start);
  record("build/empty: max_mm == 3", p.max_mm == 3);
  record("build/empty: max_bulge_dna == 2", p.max_bulge_dna == 2);
  record("build/empty: max_bulge_rna == 2", p.max_bulge_rna == 2);
}

static void test_build_empty_array_sizes() {
  auto p = make_acc(/*max_mm=*/2, /*max_bd=*/1, /*max_br=*/1).build();

  const int GL = 5, MM = 2, BD = 1, BR = 1;

  record("build/sizes: pos_mm_count.size() == guide_len",
         static_cast<int>(p.pos_mm_count.size()) == GL);
  record("build/sizes: offt_by_mm.size() == max_mm+1",
         static_cast<int>(p.offt_by_mm.size()) == MM + 1);

  record("build/sizes: ext_mm_nuc_pos.size() == max_mm+1",
         static_cast<int>(p.ext_mm_nuc_pos.size()) == MM + 1);
  record("build/sizes: ext_mm_nuc_pos[0].size() == 4",
         static_cast<int>(p.ext_mm_nuc_pos[0].size()) == 4);
  record("build/sizes: ext_mm_nuc_pos[0][0].size() == guide_len",
         static_cast<int>(p.ext_mm_nuc_pos[0][0].size()) == GL);

  record("build/sizes: ext_total_by_mm.size() == max_mm+1",
         static_cast<int>(p.ext_total_by_mm.size()) == MM + 1);
  record("build/sizes: ext_total_by_mm[0].size() == guide_len",
         static_cast<int>(p.ext_total_by_mm[0].size()) == GL);

  record("build/sizes: offt_dna.size() == max_mm+1",
         static_cast<int>(p.offt_dna.size()) == MM + 1);
  record("build/sizes: offt_dna[0].size() == max_bulge_dna",
         static_cast<int>(p.offt_dna[0].size()) == BD);

  record("build/sizes: offt_rna.size() == max_mm+1",
         static_cast<int>(p.offt_rna.size()) == MM + 1);
  record("build/sizes: offt_rna[0].size() == max_bulge_rna",
         static_cast<int>(p.offt_rna[0].size()) == BR);

  record("build/sizes: offt_complete_by_mm.size() == max_mm+1",
         static_cast<int>(p.offt_complete_by_mm.size()) == MM + 1);
  record("build/sizes: pos_mm_complete.size() == guide_len",
         static_cast<int>(p.pos_mm_complete.size()) == GL);
}

static void test_build_empty_all_zero() {
  auto p = make_acc().build();

  bool all_pos_zero = true;
  for (int v : p.pos_mm_count)
    if (v)
      all_pos_zero = false;
  for (int v : p.offt_by_mm)
    if (v)
      all_pos_zero = false;
  for (int v : p.pos_mm_complete)
    if (v)
      all_pos_zero = false;
  for (int v : p.offt_complete_by_mm)
    if (v)
      all_pos_zero = false;

  record("build/empty: all scalar counts are zero",
         p.ont_count == 0 && p.ont_count_dna == 0 && p.ont_count_rna == 0 &&
             p.ont_count_complete == 0);
  record("build/empty: all vector counts are zero", all_pos_zero);
}

// =============================================================================
// 3. MM-only push — on-target
// =============================================================================

static void test_push_mm_ontarget() {
  auto acc = make_acc();
  acc.push(ot_exact());
  auto p = acc.build();

  record("MM/ontarget: ont_count == 1", p.ont_count == 1);
  record("MM/ontarget: offt_by_mm[0] == 1", p.offt_by_mm[0] == 1);
  record("MM/ontarget: offt_by_mm[1] == 0", p.offt_by_mm[1] == 0);
  record("MM/ontarget: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
  record("MM/ontarget: ont_count_complete == 1", p.ont_count_complete == 1);
  record("MM/ontarget: offt_complete_by_mm[0] == 1",
         p.offt_complete_by_mm[0] == 1);
}

// =============================================================================
// 4. MM-only push — single mismatch at various positions
// =============================================================================

static void test_push_mm_pos0() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos0_C());
  auto p = acc.build();

  record("MM/pos0: pos_mm_count[0] == 1", p.pos_mm_count[0] == 1);
  record("MM/pos0: pos_mm_count[1] == 0", p.pos_mm_count[1] == 0);
  record("MM/pos0: offt_by_mm[1] == 1", p.offt_by_mm[1] == 1);
  record("MM/pos0: offt_by_mm[0] == 0", p.offt_by_mm[0] == 0);
  record("MM/pos0: ont_count == 0", p.ont_count == 0);
  record("MM/pos0: pos_mm_complete[0] == 1", p.pos_mm_complete[0] == 1);
  record("MM/pos0: offt_complete_by_mm[1] == 1", p.offt_complete_by_mm[1] == 1);
}

static void test_push_mm_pos2() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_A());
  auto p = acc.build();

  record("MM/pos2: pos_mm_count[2] == 1", p.pos_mm_count[2] == 1);
  record("MM/pos2: pos_mm_count[0] == 0", p.pos_mm_count[0] == 0);
  record("MM/pos2: pos_mm_count[4] == 0", p.pos_mm_count[4] == 0);
  record("MM/pos2: offt_by_mm[1] == 1", p.offt_by_mm[1] == 1);
}

static void test_push_mm_last_position() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos4_T());
  auto p = acc.build();

  record("MM/pos4: pos_mm_count[4] == 1", p.pos_mm_count[4] == 1);
  record("MM/pos4: pos_mm_count[3] == 0", p.pos_mm_count[3] == 0);
  record("MM/pos4: offt_by_mm[1] == 1", p.offt_by_mm[1] == 1);
}

static void test_push_mm_two_positions() {
  auto acc = make_acc();
  acc.push(ot_2mm_pos1G_pos3A());
  auto p = acc.build();

  record("MM/2pos: pos_mm_count[1] == 1", p.pos_mm_count[1] == 1);
  record("MM/2pos: pos_mm_count[3] == 1", p.pos_mm_count[3] == 1);
  record("MM/2pos: pos_mm_count[0] == 0", p.pos_mm_count[0] == 0);
  record("MM/2pos: pos_mm_count[2] == 0", p.pos_mm_count[2] == 0);
  record("MM/2pos: offt_by_mm[2] == 1", p.offt_by_mm[2] == 1);
  record("MM/2pos: offt_by_mm[1] == 0", p.offt_by_mm[1] == 0);
}

// =============================================================================
// 5. Nucleotide routing into ext_mm_nuc_pos
//    Nucleotide index: A=0, C=1, G=2, T=3
// =============================================================================

static void test_nucleotide_routing_A() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_nuc('A')); // guide G vs target A
  auto p = acc.build();
  record("nuc/A: ext_mm_nuc_pos[1][0][2] == 1", p.ext_mm_nuc_pos[1][0][2] == 1);
  record("nuc/A: other nucleotides at pos2 == 0",
         p.ext_mm_nuc_pos[1][1][2] == 0 && p.ext_mm_nuc_pos[1][2][2] == 0 &&
             p.ext_mm_nuc_pos[1][3][2] == 0);
  record("nuc/A: ext_total_by_mm[1][2] == 1", p.ext_total_by_mm[1][2] == 1);
}

static void test_nucleotide_routing_C() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_nuc('C'));
  auto p = acc.build();
  record("nuc/C: ext_mm_nuc_pos[1][1][2] == 1", p.ext_mm_nuc_pos[1][1][2] == 1);
  record("nuc/C: A and G slots at pos2 == 0",
         p.ext_mm_nuc_pos[1][0][2] == 0 && p.ext_mm_nuc_pos[1][2][2] == 0);
}

static void test_nucleotide_routing_G() {
  auto acc = make_acc();
  // guide[2]='G', must choose a mismatch nucleotide ≠ G, e.g. G→A then back.
  // Use pos 0 (guide A) vs G: nuc_index('G')=2.
  auto hit =
      OffTarget{"chr1", 100, Strand::Forward, "ACGTANGG", "gCGTANGG", 1, 0, 0};
  acc.push(hit);
  auto p = acc.build();
  record("nuc/G: ext_mm_nuc_pos[1][2][0] == 1", p.ext_mm_nuc_pos[1][2][0] == 1);
  record("nuc/G: other nucleotides at pos0 == 0",
         p.ext_mm_nuc_pos[1][0][0] == 0 && p.ext_mm_nuc_pos[1][1][0] == 0 &&
             p.ext_mm_nuc_pos[1][3][0] == 0);
}

static void test_nucleotide_routing_T() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_nuc('T'));
  auto p = acc.build();
  record("nuc/T: ext_mm_nuc_pos[1][3][2] == 1", p.ext_mm_nuc_pos[1][3][2] == 1);
}

// =============================================================================
// 6. Ambiguous 'N' in target — must not count as mismatch
// =============================================================================

static void test_push_N_in_target_not_mismatch() {
  auto acc = make_acc();
  acc.push(ot_N_at_pos2());
  auto p = acc.build();

  record("N/target: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
  record("N/target: ext_total_by_mm all zero",
         std::all_of(p.ext_total_by_mm[0].begin(), p.ext_total_by_mm[0].end(),
                     [](int v) { return v == 0; }));
  // The mm field in the OffTarget says 0, so it's treated as on-target.
  record("N/target: ont_count == 1", p.ont_count == 1);
  record("N/target: offt_by_mm[0] == 1", p.offt_by_mm[0] == 1);
}

static void test_push_N_in_guide_not_mismatch() {
  auto acc = make_acc();
  // Guide has 'N' at alignment pos 2 (PAM placeholder within PAM region —
  // but we test the body here). Construct explicitly with N in guide body.
  auto hit = OffTarget{"chr1",     100, Strand::Forward,
                       "ACNTANGG", // N at body pos 2
                       "ACATANGG", // A at body pos 2 — would differ from 'N'
                       0,          0,   0};
  acc.push(hit);
  auto p = acc.build();

  // toupper(N)==N, and (G != N is true, but N in guide triggers the 'gu != N'
  // guard).
  record("N/guide: pos_mm_count[2] == 0 (N in guide skips mismatch)",
         p.pos_mm_count[2] == 0);
}

// =============================================================================
// 7. DNA-bulge push
// =============================================================================

static void test_push_dna_bulge_0mm() {
  auto acc = make_acc();
  acc.push(ot_dna_bulge_pos1_0mm());
  auto p = acc.build();

  // DNA positional
  record("DNA/0mm: pos_bulge_dna[1] == 1", p.pos_bulge_dna[1] == 1);
  record("DNA/0mm: pos_bulge_dna[0] == 0", p.pos_bulge_dna[0] == 0);

  // MM-only channel untouched
  record("DNA/0mm: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
  record("DNA/0mm: offt_by_mm all zero",
         std::all_of(p.offt_by_mm.begin(), p.offt_by_mm.end(),
                     [](int v) { return v == 0; }));
  record("DNA/0mm: ont_count == 0", p.ont_count == 0);

  // DNA channel bucketing
  record("DNA/0mm: offt_dna[0][0] == 1", p.offt_dna[0][0] == 1);
  record("DNA/0mm: ont_count_dna == 1", p.ont_count_dna == 1);

  // Extended
  record("DNA/0mm: ext_dna_by_mm_pos[0][1] == 1",
         p.ext_dna_by_mm_pos[0][1] == 1);
  record("DNA/0mm: ext_dna_by_mm_pos[0][0] == 0",
         p.ext_dna_by_mm_pos[0][0] == 0);

  // Complete channel
  record("DNA/0mm: offt_complete_by_mm[0] == 1", p.offt_complete_by_mm[0] == 1);
  record("DNA/0mm: ont_count_complete == 1", p.ont_count_complete == 1);
}

static void test_push_dna_bulge_1mm() {
  auto acc = make_acc();
  acc.push(ot_dna_bulge_pos1_1mm_pos3());
  auto p = acc.build();

  // DNA positional (bulge at body pos 1)
  record("DNA/1mm: pos_bulge_dna[1] == 1", p.pos_bulge_dna[1] == 1);

  // Mismatch in DNA hit: recorded in dna channel, NOT in MM-only channel
  record("DNA/1mm: pos_mm_in_dna[3] == 1", p.pos_mm_in_dna[3] == 1);
  record("DNA/1mm: pos_mm_count[3] == 0", p.pos_mm_count[3] == 0);

  // Complete channel sees the mismatch
  record("DNA/1mm: pos_mm_complete[3] == 1", p.pos_mm_complete[3] == 1);

  // DNA bucketing
  record("DNA/1mm: offt_dna[1][0] == 1", p.offt_dna[1][0] == 1);
  record("DNA/1mm: ont_count_dna == 0", p.ont_count_dna == 0);

  // MM-only channel untouched
  record("DNA/1mm: offt_by_mm all zero",
         std::all_of(p.offt_by_mm.begin(), p.offt_by_mm.end(),
                     [](int v) { return v == 0; }));

  // Complete channel
  record("DNA/1mm: offt_complete_by_mm[1] == 1", p.offt_complete_by_mm[1] == 1);
  record("DNA/1mm: ont_count_complete == 0", p.ont_count_complete == 0);
}

// =============================================================================
// 8. RNA-bulge push
// =============================================================================

static void test_push_rna_bulge_0mm() {
  auto acc = make_acc();
  acc.push(ot_rna_bulge_pos2_0mm());
  auto p = acc.build();

  // RNA positional
  record("RNA/0mm: pos_bulge_rna[2] == 1", p.pos_bulge_rna[2] == 1);
  record("RNA/0mm: pos_bulge_rna[1] == 0", p.pos_bulge_rna[1] == 0);
  record("RNA/0mm: pos_bulge_dna all zero",
         std::all_of(p.pos_bulge_dna.begin(), p.pos_bulge_dna.end(),
                     [](int v) { return v == 0; }));

  // MM-only channel untouched
  record("RNA/0mm: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
  record("RNA/0mm: ont_count == 0", p.ont_count == 0);

  // RNA channel
  record("RNA/0mm: offt_rna[0][0] == 1", p.offt_rna[0][0] == 1);
  record("RNA/0mm: ont_count_rna == 1", p.ont_count_rna == 1);

  // Extended
  record("RNA/0mm: ext_rna_by_mm_pos[0][2] == 1",
         p.ext_rna_by_mm_pos[0][2] == 1);
  record("RNA/0mm: ext_rna_by_mm_pos[0][1] == 0",
         p.ext_rna_by_mm_pos[0][1] == 0);

  // Complete
  record("RNA/0mm: offt_complete_by_mm[0] == 1", p.offt_complete_by_mm[0] == 1);
  record("RNA/0mm: ont_count_complete == 1", p.ont_count_complete == 1);
}

static void test_push_rna_bulge_1mm() {
  auto acc = make_acc();
  acc.push(ot_rna_bulge_pos1_1mm_pos3());
  auto p = acc.build();

  // RNA positional (bulge at body pos 1)
  record("RNA/1mm: pos_bulge_rna[1] == 1", p.pos_bulge_rna[1] == 1);

  // Mismatch in RNA hit: only RNA channel, not MM-only
  record("RNA/1mm: pos_mm_in_rna[3] == 1", p.pos_mm_in_rna[3] == 1);
  record("RNA/1mm: pos_mm_count[3] == 0", p.pos_mm_count[3] == 0);
  record("RNA/1mm: pos_mm_complete[3] == 1", p.pos_mm_complete[3] == 1);

  // RNA bucketing
  record("RNA/1mm: offt_rna[1][0] == 1", p.offt_rna[1][0] == 1);
  record("RNA/1mm: ont_count_rna == 0", p.ont_count_rna == 0);

  // MM-only channel untouched
  record("RNA/1mm: offt_by_mm all zero",
         std::all_of(p.offt_by_mm.begin(), p.offt_by_mm.end(),
                     [](int v) { return v == 0; }));

  // Complete
  record("RNA/1mm: offt_complete_by_mm[1] == 1", p.offt_complete_by_mm[1] == 1);
  record("RNA/1mm: ont_count_complete == 0", p.ont_count_complete == 0);
}

// =============================================================================
// 9. PAM handling
// =============================================================================

static void test_pam_at_end_body_position_indexing() {
  // Verify that the MM at body pos 2 is recorded at index 2,
  // not shifted by the PAM.
  auto acc = make_acc(/*max_mm=*/2, 1, 1, /*pam_start=*/false);
  acc.push(ot_1mm_pos2_A());
  auto p = acc.build();

  record("PAM/end: mismatch at body pos 2, not shifted",
         p.pos_mm_count[2] == 1 && p.pos_mm_count[0] == 0 &&
             p.pos_mm_count[1] == 0);
  record("PAM/end: pam_at_start stored as false", !p.pam_at_start);
}

static void test_pam_at_start_body_position_indexing() {
  // PAM "NGG" precedes the guide body in the alignment.
  // Mismatch is at overall alignment index 5 (= pam_len=3 + body_pos=2),
  // but should be recorded at body_pos 2.
  auto acc = make_acc(/*max_mm=*/2, 1, 1, /*pam_start=*/true);
  acc.push(ot_pam_start_1mm_pos2());
  auto p = acc.build();

  record("PAM/start: mismatch recorded at body pos 2", p.pos_mm_count[2] == 1);
  record("PAM/start: body pos 0 unaffected", p.pos_mm_count[0] == 0);
  record("PAM/start: pam_at_start stored as true", p.pam_at_start);
}

static void test_pam_at_start_ontarget() {
  auto acc = make_acc(2, 1, 1, true);
  acc.push(ot_pam_start_exact());
  auto p = acc.build();
  record("PAM/start: on-target ont_count == 1", p.ont_count == 1);
  record("PAM/start: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
}

// =============================================================================
// 10. Complete channel — aggregation across all channels
// =============================================================================

static void test_complete_channel_aggregates_all() {
  auto acc = make_acc();
  acc.push(ot_exact());              // 0MM no-bulge → ont_count_complete
  acc.push(ot_dna_bulge_pos1_0mm()); // 0MM DNA → ont_count_complete
  acc.push(ot_rna_bulge_pos2_0mm()); // 0MM RNA → ont_count_complete
  acc.push(ot_1mm_pos2_A());         // 1MM no-bulge
  auto p = acc.build();

  record("complete/agg: offt_complete_by_mm[0] == 3",
         p.offt_complete_by_mm[0] == 3);
  record("complete/agg: ont_count_complete == 3 (0MM from any channel)",
         p.ont_count_complete == 3);
  record("complete/agg: offt_complete_by_mm[1] == 1",
         p.offt_complete_by_mm[1] == 1);
}

static void test_complete_pos_mm_complete_sums_channels() {
  // pos_mm_complete should sum mismatches from all channels at each position.
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_A()); // MM-only: contributes to pos_mm_count[2] and
                             // pos_mm_complete[2]
  acc.push(
      ot_dna_bulge_pos1_1mm_pos3()); // DNA+1MM: contributes to pos_mm_in_dna[3]
                                     // and pos_mm_complete[3]
  acc.push(
      ot_rna_bulge_pos1_1mm_pos3()); // RNA+1MM: contributes to pos_mm_in_rna[3]
                                     // and pos_mm_complete[3]
  auto p = acc.build();

  record("complete/pos: pos_mm_complete[2] == 1 (from MM-only hit)",
         p.pos_mm_complete[2] == 1);
  record("complete/pos: pos_mm_complete[3] == 2 (from DNA and RNA hits)",
         p.pos_mm_complete[3] == 2);
  // MM-only channel should NOT include the bulge hits' mismatch
  record("complete/pos: pos_mm_count[3] == 0", p.pos_mm_count[3] == 0);
  record("complete/pos: pos_mm_count[2] == 1", p.pos_mm_count[2] == 1);
}

// =============================================================================
// 11. Multi-hit accumulation
// =============================================================================

static void test_accumulation_multiple_mm_hits() {
  auto acc = make_acc();
  acc.push(ot_exact());           // 0MM
  acc.push(ot_1mm_pos2_A());      // 1MM @ pos2
  acc.push(ot_1mm_pos0_C());      // 1MM @ pos0
  acc.push(ot_2mm_pos1G_pos3A()); // 2MM @ pos1, pos3
  auto p = acc.build();

  record("accum: offt_by_mm[0] == 1", p.offt_by_mm[0] == 1);
  record("accum: offt_by_mm[1] == 2", p.offt_by_mm[1] == 2);
  record("accum: offt_by_mm[2] == 1", p.offt_by_mm[2] == 1);
  record("accum: pos_mm_count[0] == 1 (1MM @ pos0)", p.pos_mm_count[0] == 1);
  record("accum: pos_mm_count[1] == 1 (2MM @ pos1)", p.pos_mm_count[1] == 1);
  record("accum: pos_mm_count[2] == 1 (1MM @ pos2)", p.pos_mm_count[2] == 1);
  record("accum: pos_mm_count[3] == 1 (2MM @ pos3)", p.pos_mm_count[3] == 1);
  record("accum: pos_mm_count[4] == 0", p.pos_mm_count[4] == 0);
  record("accum: ont_count == 1", p.ont_count == 1);
}

static void test_accumulation_mixed_channels() {
  auto acc = make_acc();
  acc.push(ot_exact());                   // MM-only ONT
  acc.push(ot_dna_bulge_pos1_0mm());      // DNA ONT
  acc.push(ot_rna_bulge_pos2_0mm());      // RNA ONT
  acc.push(ot_1mm_pos0_C());              // 1MM off-target
  acc.push(ot_dna_bulge_pos1_1mm_pos3()); // DNA + 1MM
  auto p = acc.build();

  record("mixed: offt_by_mm[0] == 1 (MM-only ONT)", p.offt_by_mm[0] == 1);
  record("mixed: offt_by_mm[1] == 1 (MM-only 1MM)", p.offt_by_mm[1] == 1);
  record("mixed: offt_dna[0][0] == 1 (DNA ONT)", p.offt_dna[0][0] == 1);
  record("mixed: offt_dna[1][0] == 1 (DNA 1MM)", p.offt_dna[1][0] == 1);
  record("mixed: offt_rna[0][0] == 1 (RNA ONT)", p.offt_rna[0][0] == 1);
  record("mixed: ont_count == 1", p.ont_count == 1);
  record("mixed: ont_count_dna == 1", p.ont_count_dna == 1);
  record("mixed: ont_count_rna == 1", p.ont_count_rna == 1);
  record("mixed: ont_count_complete == 3 (all 0MM)", p.ont_count_complete == 3);
}

// =============================================================================
// 12. build() is idempotent
// =============================================================================

static void test_build_idempotent() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_A());
  acc.push(ot_dna_bulge_pos1_0mm());

  auto p1 = acc.build();
  auto p2 = acc.build(); // second call — must not alter state

  record("idempotent: offt_by_mm same after two build() calls",
         p1.offt_by_mm == p2.offt_by_mm);
  record("idempotent: pos_mm_count same", p1.pos_mm_count == p2.pos_mm_count);
  record("idempotent: offt_dna same", p1.offt_dna == p2.offt_dna);
  record("idempotent: ont_count same", p1.ont_count == p2.ont_count);
}

// =============================================================================
// 13. reset()
// =============================================================================

static void test_reset_clears_counts() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_A());
  acc.push(ot_dna_bulge_pos1_0mm());
  acc.reset();
  auto p = acc.build();

  record("reset: ont_count == 0 after reset", p.ont_count == 0);
  record("reset: ont_count_dna == 0", p.ont_count_dna == 0);
  record("reset: ont_count_complete == 0", p.ont_count_complete == 0);
  record("reset: pos_mm_count all zero",
         std::all_of(p.pos_mm_count.begin(), p.pos_mm_count.end(),
                     [](int v) { return v == 0; }));
  record("reset: offt_by_mm all zero",
         std::all_of(p.offt_by_mm.begin(), p.offt_by_mm.end(),
                     [](int v) { return v == 0; }));
  record("reset: pos_bulge_dna all zero",
         std::all_of(p.pos_bulge_dna.begin(), p.pos_bulge_dna.end(),
                     [](int v) { return v == 0; }));
  record("reset: offt_complete_by_mm all zero",
         std::all_of(p.offt_complete_by_mm.begin(), p.offt_complete_by_mm.end(),
                     [](int v) { return v == 0; }));
}

static void test_reset_preserves_geometry() {
  auto acc = make_acc(/*max_mm=*/3, /*max_bd=*/2, /*max_br=*/2);
  acc.push(ot_1mm_pos2_A());
  acc.reset();
  auto p = acc.build();

  record("reset/geometry: guide preserved", p.guide == "ACGTA");
  record("reset/geometry: guide_len preserved", p.guide_len == 5);
  record("reset/geometry: max_mm preserved", p.max_mm == 3);
  record("reset/geometry: max_bulge_dna preserved", p.max_bulge_dna == 2);
  record("reset/geometry: pos_mm_count.size() preserved",
         static_cast<int>(p.pos_mm_count.size()) == 5);
  record("reset/geometry: offt_by_mm.size() preserved",
         static_cast<int>(p.offt_by_mm.size()) == 4);
}

static void test_reset_then_push_accumulates_fresh() {
  auto acc = make_acc();
  acc.push(ot_exact());
  acc.reset();
  acc.push(ot_1mm_pos2_A());
  auto p = acc.build();

  // Should see only the post-reset push.
  record("reset/push: offt_by_mm[0] == 0 (on-target erased)",
         p.offt_by_mm[0] == 0);
  record("reset/push: offt_by_mm[1] == 1 (fresh 1MM)", p.offt_by_mm[1] == 1);
  record("reset/push: pos_mm_count[2] == 1", p.pos_mm_count[2] == 1);
}

// =============================================================================
// 14. GuideProfile value semantics
// =============================================================================

static void test_guide_profile_copy() {
  auto acc = make_acc();
  acc.push(ot_1mm_pos2_A());
  auto original = acc.build();
  auto copy = original; // copy constructor

  record("copy: offt_by_mm equal", copy.offt_by_mm == original.offt_by_mm);
  record("copy: pos_mm_count equal",
         copy.pos_mm_count == original.pos_mm_count);
  record("copy: guide equal", copy.guide == original.guide);
  record("copy: ext_mm_nuc_pos equal",
         copy.ext_mm_nuc_pos == original.ext_mm_nuc_pos);

  // Modifying copy must not affect original.
  copy.pos_mm_count[0] = 99;
  record("copy: independent from original", original.pos_mm_count[0] == 0);
}

static void test_guide_profile_move() {
  auto acc = make_acc();
  acc.push(ot_2mm_pos1G_pos3A());
  auto src = acc.build();
  const auto expected_offt = src.offt_by_mm;
  auto moved = std::move(src);

  record("move: offt_by_mm preserved", moved.offt_by_mm == expected_offt);
  record("move: max_mm preserved", moved.max_mm == 2);
  record("move: guide_len preserved", moved.guide_len == 5);
}

// =============================================================================
// 15. ext_total_by_mm consistency
//     ext_total[mm][pos] must equal Σ ext_mm_nuc_pos[mm][nuc][pos] over
//     nuc=0..3
// =============================================================================

static void test_ext_total_consistency() {
  auto acc = make_acc(/*max_mm=*/2, 1, 1, false);
  acc.push(ot_1mm_pos0_C());      // C at pos 0
  acc.push(ot_1mm_pos2_A());      // A at pos 2
  acc.push(ot_1mm_pos2_nuc('T')); // T at pos 2
  acc.push(ot_2mm_pos1G_pos3A()); // G at pos 1, A at pos 3

  auto p = acc.build();

  bool consistent = true;
  for (int mm = 0; mm <= p.max_mm; ++mm) {
    for (int pos = 0; pos < p.guide_len; ++pos) {
      int nuc_sum = 0;
      for (int nuc = 0; nuc < 4; ++nuc)
        nuc_sum += p.ext_mm_nuc_pos[mm][nuc][pos];
      if (nuc_sum != p.ext_total_by_mm[mm][pos]) {
        consistent = false;
        break;
      }
    }
  }
  record(
      "ext_total: ext_total[mm][pos] == Σ ext_nuc[mm][nuc][pos] for all mm,pos",
      consistent);
}

// =============================================================================
// 16. mm clamping: hit with mm > max_mm must not crash
// =============================================================================

static void test_mm_clamping() {
  // max_mm=1 accumulator, push a hit that reports 2 mismatches.
  auto acc = make_acc(/*max_mm=*/1, 0, 0, false);
  auto hit = OffTarget{"chr1",     100,        Strand::Forward,
                       "ACGTANGG", "aGaTANGG", // 2 mismatches
                       2,          0,          0};
  bool threw = false;
  try {
    acc.push(hit);
  } catch (...) {
    threw = true;
  }
  record("clamp: push with mm > max_mm does not throw", !threw);

  auto p = acc.build();
  // mm_idx clamped to 1; offt_by_mm[1] should be 1 (not out-of-range).
  record("clamp: hit bucketed at clamped mm_idx", p.offt_by_mm[1] == 1);
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_profile_data ===\n\n";

  std::cout << "-- Construction --\n";
  test_construction_valid();
  test_construction_invalid_guide_len_zero();
  test_construction_invalid_guide_len_negative();
  test_construction_invalid_pam_len_negative();
  test_construction_invalid_max_mm_negative();
  test_construction_invalid_max_bulge_dna_negative();
  test_construction_invalid_max_bulge_rna_negative();
  test_construction_boundary_zero_mm_no_bulge();

  std::cout << "\n-- build() on empty accumulator --\n";
  test_build_empty_geometry();
  test_build_empty_array_sizes();
  test_build_empty_all_zero();

  std::cout << "\n-- MM-only push --\n";
  test_push_mm_ontarget();
  test_push_mm_pos0();
  test_push_mm_pos2();
  test_push_mm_last_position();
  test_push_mm_two_positions();

  std::cout << "\n-- Nucleotide routing --\n";
  test_nucleotide_routing_A();
  test_nucleotide_routing_C();
  test_nucleotide_routing_G();
  test_nucleotide_routing_T();

  std::cout << "\n-- Ambiguous N --\n";
  test_push_N_in_target_not_mismatch();
  test_push_N_in_guide_not_mismatch();

  std::cout << "\n-- DNA-bulge push --\n";
  test_push_dna_bulge_0mm();
  test_push_dna_bulge_1mm();

  std::cout << "\n-- RNA-bulge push --\n";
  test_push_rna_bulge_0mm();
  test_push_rna_bulge_1mm();

  std::cout << "\n-- PAM handling --\n";
  test_pam_at_end_body_position_indexing();
  test_pam_at_start_body_position_indexing();
  test_pam_at_start_ontarget();

  std::cout << "\n-- Complete channel --\n";
  test_complete_channel_aggregates_all();
  test_complete_pos_mm_complete_sums_channels();

  std::cout << "\n-- Multi-hit accumulation --\n";
  test_accumulation_multiple_mm_hits();
  test_accumulation_mixed_channels();

  std::cout << "\n-- build() idempotency --\n";
  test_build_idempotent();

  std::cout << "\n-- reset() --\n";
  test_reset_clears_counts();
  test_reset_preserves_geometry();
  test_reset_then_push_accumulates_fresh();

  std::cout << "\n-- GuideProfile value semantics --\n";
  test_guide_profile_copy();
  test_guide_profile_move();

  std::cout << "\n-- ext_total consistency --\n";
  test_ext_total_consistency();

  std::cout << "\n-- mm clamping --\n";
  test_mm_clamping();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}