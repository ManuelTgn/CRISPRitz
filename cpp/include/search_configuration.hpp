#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace crispritz {

// =========================================================================
// OutputFormat
// =========================================================================

/**
 * @brief Serialization format for off-target search results.
 *
 * Scoped enum so an invalid format is unrepresentable at compile time,
 * mirroring the design of @c Strand in off_target.hpp.  The underlying
 * values are stable and may be used as array indices, but must not be
 * reordered if they are ever persisted.
 */
enum class OutputFormat : uint8_t {
  Tsv = 0,     ///< Tab-separated values (one row per off-target).
  Targets = 1, ///< Legacy CRISPRitz "targets" column layout.
};

/**
 * @brief Return the lowercase canonical name of an @c OutputFormat.
 *
 * @param fmt  Output format value.
 * @return     "tsv" or "targets" — a view into static storage.
 * @complexity O(1).
 */
[[nodiscard]] std::string_view to_string(OutputFormat fmt) noexcept;

/**
 * @brief Parse a format name (case-sensitive, lowercase) into the enum.
 *
 * @param name  "tsv" or "targets".
 * @return      The corresponding @c OutputFormat.
 * @throws std::invalid_argument if @p name is not a recognised format.
 * @complexity O(1).
 */
[[nodiscard]] OutputFormat output_format_from_string(std::string_view name);

// =========================================================================
// OutputMode
// =========================================================================

/**
 * @brief Controls which output files are produced after a search run.
 *
 * Mirrors the legacy @c -t flag behaviour in @c searchOnTST.cpp:
 *   @c "r"  → TargetsOnly  (only @c .targets.txt / @c .tsv)
 *   @c "p"  → ProfileOnly  (only @c .profile.xls family)
 *   default → Both
 *
 * The mode is resolved once at the CLI / Python entry point and stored in
 * @c SearchConfiguration so that the orchestration layer can branch on it
 * via the convenience helpers @c write_targets() / @c write_profile()
 * without inspecting raw strings deep in the pipeline.
 *
 * ## Effect on the accumulation loop
 * When @c TargetsOnly, no @c ProfileAccumulator objects are constructed
 * and @c ProfileAccumulator::push() is never called — the profiling
 * infrastructure is entirely skipped.  When @c ProfileOnly, the
 * @c OutputWriter::Session is not opened.  When @c Both, both paths
 * are active in the same per-hit loop.
 */
enum class OutputMode : uint8_t {
  TargetsOnly = 0, ///< Write targets table only.
  ProfileOnly = 1, ///< Write profile files only (.profile.xls, …).
  Both = 2,        ///< Write both targets and all profile files (default).
};

/**
 * @brief Return the lowercase canonical name of an @c OutputMode.
 *
 * @param mode  Output mode value.
 * @return      "targets", "profile", or "both" — a view into static storage.
 * @complexity  O(1).
 */
[[nodiscard]] std::string_view to_string(OutputMode mode) noexcept;

/**
 * @brief Parse a mode name (case-sensitive, lowercase) into the enum.
 *
 * @param name  "targets", "profile", or "both".
 * @return      The corresponding @c OutputMode.
 * @throws std::invalid_argument if @p name is not recognised.
 * @complexity  O(1).
 */
[[nodiscard]] OutputMode output_mode_from_string(std::string_view name);

// =========================================================================
// SearchConfiguration
// =========================================================================

/**
 * @brief Immutable bundle of the search-time parameters a user controls.
 *
 * One @c SearchConfiguration is built once (from CLI arguments or the
 * Python entry point), validated, and passed by @c const& into the search
 * engine.  Because it is immutable and carries no mutable state, a single
 * instance can be shared across all worker threads without synchronisation.
 *
 * ## Scope: search-time user decisions and output mode
 * This object deliberately excludes the PAM specification and guide length.
 * Those are *indexing-time* properties: the PAM was already applied during
 * the genome-index phase (only sites presenting a valid PAM were extracted
 * into the TST), and the resulting geometry — guide length, PAM limit, PAM
 * orientation — is persisted in each @c .bin partition header.  The search
 * layer obtains that geometry from the loaded index, not from the user, so
 * replicating it here would create two sources of truth that can disagree.
 *
 * ## Stored parameters
 *   - @c max_mismatches — maximum substitution mismatches (>= 0).
 *   - @c max_bulges_dna — maximum DNA-bulge bases (>= 0).
 *   - @c max_bulges_rna — maximum RNA-bulge bases (>= 0).
 *   - @c threads        — worker thread count (>= 1).
 *   - @c output_format  — result serialization format.
 *   - @c output_mode    — which output files to produce.
 *
 * ## Invariants enforced here (standalone-checkable)
 *   - max_mismatches >= 0
 *   - max_bulges_dna >= 0
 *   - max_bulges_rna >= 0
 *   - threads        >= 1
 *
 * ## Invariant NOT enforced here (deferred — see note)
 * The relationship @c max_total_edits() <= guide_length is an
 * index-dependent constraint: an alignment cannot contain more edits than
 * the guide has positions.  Because @c guide_length is read from the index
 * header at load time, this check belongs in the index loader / search
 * engine once the header is available, not in this standalone object.
 * @c max_total_edits() is exposed precisely so that layer can perform the
 * comparison with a single call.
 */
class SearchConfiguration {
public:
  /**
   * @brief Construct and validate the standalone-checkable parameters.
   *
   * @param max_mismatches Maximum substitution mismatches (>= 0).
   * @param max_bulges_dna Maximum DNA-bulge bases (>= 0).
   * @param max_bulges_rna Maximum RNA-bulge bases (>= 0).
   * @param threads        Worker thread count (>= 1).
   * @param output_format  Result serialization format
   *                       (default: OutputFormat::Tsv).
   * @param output_mode    Which output files to produce
   *                       (default: OutputMode::Both).
   *
   * @throws std::invalid_argument if any numeric invariant is violated.
   *
   * @note This constructor does NOT validate the edit budget against the
   *       guide length; that check is deferred to the index loader (see
   *       class documentation).
   */
  SearchConfiguration(int max_mismatches, int max_bulges_dna,
                      int max_bulges_rna, int threads,
                      OutputFormat output_format = OutputFormat::Tsv,
                      OutputMode output_mode = OutputMode::Both);

  // ---- Accessors (all noexcept) ----------------------------------------

  /** @return Maximum substitution mismatches. */
  [[nodiscard]] int max_mismatches() const noexcept { return max_mismatches_; }

  /** @return Maximum DNA-bulge bases. */
  [[nodiscard]] int max_bulges_dna() const noexcept { return max_bulges_dna_; }

  /** @return Maximum RNA-bulge bases. */
  [[nodiscard]] int max_bulges_rna() const noexcept { return max_bulges_rna_; }

  /** @return Worker thread count (>= 1). */
  [[nodiscard]] int threads() const noexcept { return threads_; }

  /** @return Result serialization format. */
  [[nodiscard]] OutputFormat output_format() const noexcept {
    return output_format_;
  }

  /** @return Which output files to produce. */
  [[nodiscard]] OutputMode output_mode() const noexcept { return output_mode_; }

  // ---- Derived properties -----------------------------------------------

  /**
   * @brief Total bulge budget: DNA bulges + RNA bulges.
   * @return max_bulges_dna() + max_bulges_rna().
   */
  [[nodiscard]] int max_bulges_total() const noexcept {
    return max_bulges_dna_ + max_bulges_rna_;
  }

  /**
   * @brief Total edit budget the search recursion may spend.
   *
   * Exposed so the index loader can validate it against the index's
   * guide length once the header is available (see class documentation).
   *
   * @return max_mismatches() + max_bulges_total().
   */
  [[nodiscard]] int max_total_edits() const noexcept {
    return max_mismatches_ + max_bulges_total();
  }

  /**
   * @brief True when no bulges are allowed (DNA and RNA budgets both 0).
   *
   * The search engine can use this to select the cheaper mismatch-only
   * traversal kernel.
   *
   * @return max_bulges_total() == 0.
   */
  [[nodiscard]] bool bulges_disabled() const noexcept {
    return max_bulges_total() == 0;
  }

  // ---- Output-mode convenience helpers ---------------------------------
  //
  // These two predicates are the single branch point for the orchestration
  // loop.  All code that needs to decide "do I open a targets Session?"
  // or "do I create ProfileAccumulators?" calls one of these rather than
  // comparing output_mode() directly, so the enum values stay encapsulated
  // here.

  /**
   * @brief True when the targets table should be written.
   * @return @c true for TargetsOnly and Both.
   */
  [[nodiscard]] bool write_targets() const noexcept {
    return output_mode_ == OutputMode::TargetsOnly ||
           output_mode_ == OutputMode::Both;
  }

  /**
   * @brief True when profile files should be written.
   * @return @c true for ProfileOnly and Both.
   */
  [[nodiscard]] bool write_profile() const noexcept {
    return output_mode_ == OutputMode::ProfileOnly ||
           output_mode_ == OutputMode::Both;
  }

private:
  int max_mismatches_;
  int max_bulges_dna_;
  int max_bulges_rna_;
  int threads_;
  OutputFormat output_format_;
  OutputMode output_mode_;
};

} // namespace crispritz