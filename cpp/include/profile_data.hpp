#pragma once

/**
 * @file profile_data.hpp
 * @brief Per-guide profiling statistics and the accumulator that builds them.
 *
 * ## Design overview
 *
 * The legacy CRISPRitz code maintained several large multi-dimensional matrices
 * (profiling, ext_profiling, profiling_dna_mm, …) that were updated per-hit
 * *inside* the search loop, across multiple OpenMP threads.  In the refactored
 * architecture those matrices are replaced by two types defined here:
 *
 *   - @c GuideProfile  — a plain-data value type that owns all profiling
 *     statistics for one guide after aggregation.  One instance corresponds to
 *     one row (or block) in each of the five legacy profile files.
 *
 *   - @c ProfileAccumulator — a stateful, single-guide accumulator that ingests
 *     @c OffTarget objects one at a time and builds the corresponding
 *     @c GuideProfile counters.  It replaces the @c detailedOutputFast* family
 *     from @c detailedOutput.cpp.
 *
 * ## Accumulation vs. search separation
 *
 * @c ProfileAccumulator sits entirely *outside* the search layer.  It consumes
 * @c OffTarget objects that have already been produced by @c TSTSearcher — it
 * never reaches into the TST or the alignment kernel.  This means:
 *
 *   - @c TSTSearcher is unmodified and profile-unaware.
 *   - Accumulation and the targets-table flush session share the same per-hit
 *     loop in the orchestration layer (Python or search_runner.cpp) and can
 *     therefore be composed cheaply without a second pass over the hits.
 *   - Memory cost is dominated by the accumulated counters (~3–4 KB per guide)
 *     not by the full hit list, so all guides' accumulators can live in memory
 *     simultaneously even for large guide sets.
 *
 * ## Gap conventions in OffTarget alignment strings
 *
 *   - @c grna()[i] == '-'   → DNA bulge at alignment column @c i
 *                             (target has an extra base; guide skips it).
 *   - @c target()[i] == '-' → RNA bulge at alignment column @c i
 *                             (guide has an extra base; target skips it).
 *   - Both non-gap, uppercase-equal → match.
 *   - Both non-gap, uppercase-differ, neither 'N' → substitution mismatch.
 *
 * ## Output files produced from GuideProfile
 *
 * | File                     | Fields used                                    |
 * |--------------------------|------------------------------------------------|
 * | @c .profile.xls          | pos_mm_count, ont_count, offt_by_mm            |
 * | @c .extended_profile.xls | ext_mm_nuc_pos, ext_total_by_mm,               |
 * |                          | ext_dna_by_mm_pos, ext_rna_by_mm_pos           |
 * | @c .profile_dna.xls      | pos_bulge_dna, pos_mm_in_dna, offt_dna,        |
 * |                          | ont_count_dna                                  |
 * | @c .profile_rna.xls      | pos_bulge_rna, pos_mm_in_rna, offt_rna,        |
 * |                          | ont_count_rna                                  |
 * | @c .profile_complete.xls | pos_mm_complete, ont_count_complete, | | |
 * offt_complete_by_mm                            |
 */

#include "offtarget.hpp"

#include <string>
#include <vector>

namespace crispritz {

// =============================================================================
// GuideProfile
// =============================================================================

/**
 * @brief All profiling statistics for a single guide, post-accumulation.
 *
 * A @c GuideProfile is produced by @c ProfileAccumulator::build() after all
 * @c OffTarget hits for a guide have been pushed.  It is a pure-data value
 * type: all fields are public, there are no invariants to maintain, and it
 * participates in regular C++ value semantics (copy, move).
 *
 * ## Indexing conventions
 *
 *   - Position indices (0-based) refer to the *guide body* — PAM bases are
 *     not counted, even though the PAM appears in the alignment string.
 *     The accumulator strips PAM columns before updating position arrays.
 *
 *   - Mismatch-count indices run 0 … max_mm (inclusive).  Index 0 denotes
 *     on-target hits (zero mismatches, zero bulge); indices ≥ 1 denote
 *     off-target hits with exactly that many mismatches.
 *
 *   - Bulge-size indices run 0 … max_bulge-1, i.e. bulge size 1 → index 0,
 *     bulge size k → index k-1.  This matches the legacy matrix layout.
 *
 *   - Nucleotide indices: A=0, C=1, G=2, T=3 (the off-target base at the
 *     mismatch position, matching the legacy @c matrixprofiling convention).
 */
struct GuideProfile {
  // ------------------------------------------------------------------
  // Identity / geometry (set by ProfileAccumulator constructor)
  // ------------------------------------------------------------------

  /// Guide body sequence, PAM bases excluded.
  std::string guide;

  /// Length of the guide body (== guide.size(), stored for convenience).
  int guide_len{0};

  /// Number of PAM bases appended/prepended (needed for header row Ns).
  int pam_len{0};

  /// @c true when the PAM precedes the guide body in the alignment string.
  bool pam_at_start{false};

  /// Maximum mismatches the search was configured with.
  int max_mm{0};

  /// Maximum DNA-bulge bases the search was configured with.
  int max_bulge_dna{0};

  /// Maximum RNA-bulge bases the search was configured with.
  int max_bulge_rna{0};

  // ------------------------------------------------------------------
  // File 1: .profile.xls  (mismatch-only channel)
  // ------------------------------------------------------------------

  /**
   * @brief Per-position mismatch event count.
   *
   * @c pos_mm_count[i] is the number of off-target hits (across all
   * channels) in which alignment column @c i carried a substitution
   * mismatch.  Size: @c guide_len.
   */
  std::vector<int> pos_mm_count;

  /**
   * @brief On-target hit count (MM-only channel, 0 mismatches, 0 bulge).
   *
   * Corresponds to the @c ONT column of the legacy @c .profile.xls file.
   */
  int ont_count{0};

  /**
   * @brief Per-mismatch-count off-target hit tally (MM-only, no bulge).
   *
   * @c offt_by_mm[n] is the number of hits with exactly @c n mismatches
   * and zero bulge bases.  Index 0 equals @c ont_count.
   * Size: @c max_mm + 1.
   */
  std::vector<int> offt_by_mm;

  // ------------------------------------------------------------------
  // File 2: .extended_profile.xls
  // ------------------------------------------------------------------

  /**
   * @brief Extended per-(mm, nucleotide, position) hit count.
   *
   * @c ext_mm_nuc_pos[mm][nuc][pos] is the number of hits that had
   * exactly @c mm total mismatches and carried nucleotide @c nuc
   * (A=0, C=1, G=2, T=3) at guide body position @c pos in the
   * off-target sequence.
   * Dimensions: (max_mm+1) × 4 × guide_len.
   */
  std::vector<std::vector<std::vector<int>>> ext_mm_nuc_pos;

  /**
   * @brief Total mismatch events per (mm threshold, guide position).
   *
   * @c ext_total_by_mm[mm][pos] = sum over nucleotides of
   * @c ext_mm_nuc_pos[mm][nuc][pos].
   * Dimensions: (max_mm+1) × guide_len.
   */
  std::vector<std::vector<int>> ext_total_by_mm;

  /**
   * @brief DNA-bulge events per (mm threshold, guide body position).
   *
   * @c ext_dna_by_mm_pos[mm][pos] is the number of DNA-bulge events
   * at body position @c pos among hits that had exactly @c mm
   * mismatches.
   * Dimensions: (max_mm+1) × guide_len.
   */
  std::vector<std::vector<int>> ext_dna_by_mm_pos;

  /**
   * @brief RNA-bulge events per (mm threshold, guide body position).
   *
   * Same structure as @c ext_dna_by_mm_pos but for RNA bulges.
   * Dimensions: (max_mm+1) × guide_len.
   */
  std::vector<std::vector<int>> ext_rna_by_mm_pos;

  // ------------------------------------------------------------------
  // File 3: .profile_dna.xls  (DNA-bulge channel)
  // ------------------------------------------------------------------

  /**
   * @brief DNA-bulge event count per guide body position.
   *
   * @c pos_bulge_dna[i] is the total number of DNA-bulge events
   * observed at body position @c i across all hits.  Size: @c guide_len.
   */
  std::vector<int> pos_bulge_dna;

  /**
   * @brief Mismatch event count per body position in DNA-bulge alignments.
   *
   * @c pos_mm_in_dna[i] counts mismatches at body position @c i that
   * occurred within a DNA-bulge alignment (i.e. @c bulge_dna() > 0).
   * Size: @c guide_len.
   */
  std::vector<int> pos_mm_in_dna;

  /**
   * @brief DNA-bulge off-target count per (mm count, bulge size).
   *
   * @c offt_dna[mm][b] is the number of hits with exactly @c mm
   * mismatches and @c b+1 DNA-bulge bases.
   * Dimensions: (max_mm+1) × max_bulge_dna.
   */
  std::vector<std::vector<int>> offt_dna;

  /**
   * @brief On-target count for the DNA-bulge channel.
   *
   * Hits with 0 mismatches and exactly 1 or more DNA-bulge bases.
   */
  int ont_count_dna{0};

  // ------------------------------------------------------------------
  // File 4: .profile_rna.xls  (RNA-bulge channel)
  // ------------------------------------------------------------------

  /** @brief RNA-bulge event count per guide body position. Size: guide_len. */
  std::vector<int> pos_bulge_rna;

  /**
   * @brief Mismatch count per body position in RNA-bulge alignments.
   * Size: @c guide_len.
   */
  std::vector<int> pos_mm_in_rna;

  /**
   * @brief RNA-bulge off-target count per (mm count, bulge size).
   * Dimensions: (max_mm+1) × max_bulge_rna.
   */
  std::vector<std::vector<int>> offt_rna;

  /** @brief On-target count for the RNA-bulge channel (0 mm, ≥1 RNA bulge). */
  int ont_count_rna{0};

  // ------------------------------------------------------------------
  // File 5: .profile_complete.xls  (all-channel combined)
  // ------------------------------------------------------------------

  /**
   * @brief Combined off-target count per total-mismatch count.
   *
   * @c offt_complete_by_mm[n] is the number of hits with exactly @c n
   * mismatches across all bulge channels (MM-only + DNA + RNA).
   * Size: @c max_mm + 1.
   */
  std::vector<int> offt_complete_by_mm;

  /**
   * @brief On-target count across all channels (0 mismatches, any bulge).
   */
  int ont_count_complete{0};

  /**
   * @brief All-channel mismatch event count per guide body position.
   *
   * Sum of mismatch events at each body position across MM-only, DNA-bulge,
   * and RNA-bulge hits.  Size: @c guide_len.
   */
  std::vector<int> pos_mm_complete;
};

// =============================================================================
// ProfileAccumulator
// =============================================================================

/**
 * @brief Stateful per-guide accumulator that builds a @c GuideProfile.
 *
 * ## Lifecycle
 * 1. Construct once per guide with the guide's sequence and search geometry.
 * 2. Call @c push() for every @c OffTarget hit belonging to that guide,
 *    in any order.
 * 3. Call @c build() to obtain the finished @c GuideProfile (or @c reset()
 *    to clear counters and reuse the accumulator for a fresh search).
 *
 * ## Thread safety
 * A single @c ProfileAccumulator instance is **not** thread-safe.  In the
 * partition-parallel design each partition is processed by one Python-thread
 * at a time; hits from that partition are pushed sequentially, then the
 * guide's merged accumulator is updated in a critical section in Python.
 * No locking is needed inside this class.
 *
 * @complexity
 *   push(): O(guide_len + pam_len) per call (one pass over the alignment).
 *   build(): O(guide_len × max_mm) (copying and summing counters).
 *   reset(): O(guide_len × max_mm) (zero-filling).
 *   Memory: O(guide_len × max_mm) integers.
 */
class ProfileAccumulator {
public:
  /**
   * @brief Construct the accumulator for one guide.
   *
   * All internal counter arrays are zero-initialised.
   *
   * @param guide        Guide body sequence (PAM bases excluded).
   * @param guide_len    Length of the guide body (== guide.size()).
   * @param pam_len      Number of PAM bases (appended or prepended).
   * @param max_mm       Maximum mismatches the search was configured with.
   * @param max_bulge_dna Maximum DNA-bulge bases.
   * @param max_bulge_rna Maximum RNA-bulge bases.
   * @param pam_at_start @c true when the PAM precedes the guide body in
   *                     the alignment strings (e.g. SpCas9 with 5′ PAM).
   *
   * @throws std::invalid_argument if guide_len <= 0, pam_len < 0,
   *         max_mm < 0, max_bulge_dna < 0, or max_bulge_rna < 0.
   */
  ProfileAccumulator(std::string guide, int guide_len, int pam_len, int max_mm,
                     int max_bulge_dna, int max_bulge_rna, bool pam_at_start);

  /**
   * @brief Ingest one off-target hit and update all counters.
   *
   * Walks @p ot.grna() and @p ot.target() character by character.
   * PAM columns are identified by position (leading or trailing, according
   * to @c pam_at_start_) and skipped for body-position indexing.
   * Per alignment column:
   *   - @c grna[i]=='-' : DNA bulge — increments DNA counters.
   *   - @c target[i]=='-': RNA bulge — increments RNA counters.
   *   - Both non-gap, different (non-N): mismatch — increments
   *     per-position and per-nucleotide counters.
   *
   * At the end of the alignment, the hit is bucketed into the appropriate
   * @c offt_by_mm / @c offt_dna / @c offt_rna cell.
   *
   * @param ot  The off-target hit to accumulate.
   * @complexity O(guide_len + pam_len).
   */
  void push(const OffTarget &ot);

  /**
   * @brief Produce the finished @c GuideProfile from the accumulated counters.
   *
   * Does not reset the accumulator; call @c reset() explicitly if the
   * accumulator is to be reused.
   *
   * @return A freshly constructed @c GuideProfile by value.
   * @complexity O(guide_len × max_mm).
   */
  [[nodiscard]] GuideProfile build() const;

  /**
   * @brief Zero all counters, preserving the guide identity and geometry.
   *
   * Allows the accumulator to be reused for a new search over the same
   * guide without reallocation.
   *
   * @complexity O(guide_len × max_mm).
   */
  void reset();

private:
  // ---- Identity / geometry -------------------------------------------
  std::string guide_;
  int guide_len_;
  int pam_len_;
  int max_mm_;
  int max_bulge_dna_;
  int max_bulge_rna_;
  bool pam_at_start_;

  // ---- File 1 counters -----------------------------------------------
  std::vector<int> pos_mm_count_; // [guide_len]
  int ont_count_{0};
  std::vector<int> offt_by_mm_; // [max_mm+1]

  // ---- File 2 counters -----------------------------------------------
  // [max_mm+1][4][guide_len]
  std::vector<std::vector<std::vector<int>>> ext_mm_nuc_pos_;
  // [max_mm+1][guide_len]
  std::vector<std::vector<int>> ext_total_by_mm_;
  std::vector<std::vector<int>> ext_dna_by_mm_pos_;
  std::vector<std::vector<int>> ext_rna_by_mm_pos_;

  // ---- File 3 counters -----------------------------------------------
  std::vector<int> pos_bulge_dna_;         // [guide_len]
  std::vector<int> pos_mm_in_dna_;         // [guide_len]
  std::vector<std::vector<int>> offt_dna_; // [max_mm+1][max_bulge_dna]
  int ont_count_dna_{0};

  // ---- File 4 counters -----------------------------------------------
  std::vector<int> pos_bulge_rna_;
  std::vector<int> pos_mm_in_rna_;
  std::vector<std::vector<int>> offt_rna_; // [max_mm+1][max_bulge_rna]
  int ont_count_rna_{0};

  // ---- File 5 counters -----------------------------------------------
  std::vector<int> offt_complete_by_mm_; // [max_mm+1]
  int ont_count_complete_{0};
  std::vector<int> pos_mm_complete_; // [guide_len]

  // ---- Internal helpers ----------------------------------------------

  /**
   * @brief Map a nucleotide character to its index (A=0, C=1, G=2, T=3).
   *
   * The character is converted to uppercase before mapping.  'N' and any
   * other non-ACGT character return -1 so callers can skip ambiguous bases.
   *
   * @param c  Nucleotide character.
   * @return   0–3 for A/C/G/T; -1 otherwise.
   */
  [[nodiscard]] static int nuc_index(char c) noexcept;
};

} // namespace crispritz