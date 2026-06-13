#pragma once

#include <cstdint>
#include <functional> // std::hash
#include <string>
#include <vector>

namespace crispritz {

// =========================================================================
// Strand
// =========================================================================

/**
 * @brief Strand orientation of a genomic off-target site.
 *
 * The underlying type is @c char so that a lossless cast to the
 * single-character TSV representation is always available via
 * @c to_char().  The enum is scoped to prevent implicit conversion
 * from a raw @c char at call sites.
 */
enum class Strand : char { Forward = '+', Reverse = '-' };

/**
 * @brief Convert a Strand value to its single-character representation.
 *
 * @param s  Strand enum value.
 * @return   '+' for Forward, '-' for Reverse.
 */
constexpr char to_char(Strand s) noexcept { return static_cast<char>(s); }

/**
 * @brief Parse a strand character into a @c Strand enum.
 *
 * @param c  Must be '+' or '-'.
 * @return   The corresponding Strand value.
 * @throws   std::invalid_argument if @p c is any other character.
 */
Strand strand_from_char(char c);

// =========================================================================
// OffTarget
// =========================================================================

/**
 * @brief Represents a single CRISPR off-target site produced by TST search.
 *
 * An @c OffTarget captures the fully resolved genomic location, aligned
 * sequences, and edit-distance breakdown for one candidate off-target locus.
 * It is the canonical in-memory domain object produced by @c TSTSearcher
 * and consumed by the serialization and Python scoring layers.
 *
 * ## What this class does NOT store
 * - CFD score or any other post-search score.  These are computed in the
 *   Python layer and appended to the output file separately.
 * - PAM sequence as a separate field.  The PAM bases appear at the end of
 *   @p target (or the beginning, for upstream PAMs); callers that need the
 *   PAM slice can extract it from @c target() using @c pam_limit.
 *
 * ## Identity semantics
 * Two @c OffTarget objects are considered equal when they represent the
 * same genomic window: same chromosome, same 1-based position, same strand,
 * and same aligned target sequence.  The guide RNA sequence and the numeric
 * edit-distance breakdown are intentionally excluded from equality so that
 * deduplication across index partitions works correctly even when the same
 * site is reached by different TST traversal paths.
 *
 * ## Ordering semantics
 * @c operator< sorts by chromosome lexicographically, then position, then
 * strand character, then total edit distance.  For natural chromosome order
 * (chr1 < chr2 < chr10) callers should supply a custom comparator at the
 * output layer; lexicographic order places chr10 before chr2.
 *
 * ## Note on field naming vs. Python API
 * The Python @c OffTarget class calls the genomic target field @e spacer
 * for historical reasons.  This C++ class uses the more accurate term
 * @e target (the genomic sequence the guide RNA aligns to).  The Python
 * integration layer maps @c target() → @c spacer in the final TSV.
 */
class OffTarget {
public:
  // =====================================================================
  // Construction
  // =====================================================================

  /**
   * @brief Construct an @c OffTarget from all resolved search fields.
   *
   * String arguments are accepted by value and @em moved into the
   * object; pass @c std::move(s) for temporaries and named locals
   * that will not be used again.
   *
   * @param chrom      Chromosome / contig name (e.g. @c "chr1").
   *                   Must be non-empty.
   * @param pos        1-based genomic start position.  Must be > 0.
   * @param strand     Strand orientation.
   * @param grna       Aligned guide RNA sequence.  Includes PAM-placeholder
   *                   @c N characters at the PAM positions.  Mismatch
   *                   positions are uppercase; @c '-' marks bulge gaps.
   *                   Must be non-empty.
   * @param target     Aligned genomic target sequence.  Includes actual PAM
   *                   bases.  Mismatch positions are lowercase (legacy
   *                   convention from the C++ search); @c '-' marks bulge
   *                   gaps.  Must be non-empty.
   * @param mismatches Number of substitution mismatches.  Must be >= 0.
   * @param bulge_dna  Number of DNA-bulge bases (gap in guide, extra base
   *                   in target).  Must be >= 0.
   * @param bulge_rna  Number of RNA-bulge bases (gap in target, extra base
   *                   in guide).  Must be >= 0.
   *
   * @throws std::invalid_argument if any of the above constraints is
   *         violated.
   */
  OffTarget(std::string chrom, int32_t pos, Strand strand, std::string grna,
            std::string target, int mismatches, int bulge_dna, int bulge_rna);

  // =====================================================================
  // Accessors  (all noexcept — returning references to stored members)
  // =====================================================================

  /** @return Chromosome / contig name. */
  [[nodiscard]] const std::string &chrom() const noexcept { return chrom_; }

  /** @return 1-based genomic start position. */
  [[nodiscard]] int32_t pos() const noexcept { return pos_; }

  /** @return Strand orientation (Forward '+' or Reverse '-'). */
  [[nodiscard]] Strand strand() const noexcept { return strand_; }

  /**
   * @return Aligned guide RNA sequence with PAM-placeholder Ns.
   *         Mismatch positions uppercase; bulge gaps marked with '-'.
   */
  [[nodiscard]] const std::string &grna() const noexcept { return grna_; }

  /**
   * @return Aligned genomic target sequence including PAM bases.
   *         Mismatch positions lowercase; bulge gaps marked with '-'.
   */
  [[nodiscard]] const std::string &target() const noexcept { return target_; }

  /** @return Number of substitution mismatches. */
  [[nodiscard]] int mismatches() const noexcept { return mismatches_; }

  /** @return Number of DNA-bulge bases (gap in guide, extra in target). */
  [[nodiscard]] int bulge_dna() const noexcept { return bulge_dna_; }

  /** @return Number of RNA-bulge bases (gap in target, extra in guide). */
  [[nodiscard]] int bulge_rna() const noexcept { return bulge_rna_; }

  // =====================================================================
  // Derived properties
  // =====================================================================

  /**
   * @brief Total edit distance: mismatches + DNA-bulge bases + RNA-bulge
   *        bases.
   *
   * This is the primary filter criterion during search and the secondary
   * sort key in output files.
   *
   * @return mismatches() + bulge_dna() + bulge_rna().
   */
  [[nodiscard]] int total_edit_distance() const noexcept {
    return mismatches_ + bulge_dna_ + bulge_rna_;
  }

  /**
   * @brief True when the alignment contains at least one bulge position.
   * @return @c true if bulge_dna() > 0 or bulge_rna() > 0.
   */
  [[nodiscard]] bool has_bulge() const noexcept {
    return bulge_dna_ > 0 || bulge_rna_ > 0;
  }

  /**
   * @brief Human-readable bulge classification, derived from the bulge
   *        counts.
   *
   * @return @c "X"       — no bulge (bulge_dna == 0 && bulge_rna == 0).
   *         @c "DNA"     — DNA bulge only.
   *         @c "RNA"     — RNA bulge only.
   *         @c "DNA,RNA" — both types present.
   */
  [[nodiscard]] std::string bulge_type() const;

  // =====================================================================
  // Comparison operators
  // =====================================================================

  /**
   * @brief Locus equality: same chrom, pos, strand, and target sequence.
   *
   * The guide RNA sequence and numeric edit-distance breakdown are
   * intentionally excluded.  This allows deduplication of the same
   * genomic site found across multiple index partitions or via different
   * traversal paths, regardless of how the alignment is represented.
   */
  [[nodiscard]] bool operator==(const OffTarget &other) const noexcept;

  /** @return !(*this == other) */
  [[nodiscard]] bool operator!=(const OffTarget &other) const noexcept;

  /**
   * @brief Genomic coordinate ordering.
   *
   * Sort key in priority order:
   *   1. @c chrom  — lexicographic (not natural; chr10 < chr2).
   *   2. @c pos    — ascending.
   *   3. @c strand — '+' (43) < '-' (45) by ASCII value.
   *   4. @c total_edit_distance() — ascending tiebreaker.
   *
   * This ordering is a strict weak order and is safe for @c std::sort
   * and @c std::set.  For natural chromosome order, supply a custom
   * comparator at the output layer.
   */
  [[nodiscard]] bool operator<(const OffTarget &other) const noexcept;

  // =====================================================================
  // Serialization
  // =====================================================================

  /**
   * @brief Produce a separated-value row for the result file.
   *
   * Column order (matches @c tsv_header()):
   *   chrom, pos, strand, grna, target, mismatches,
   *   bulge_dna, bulge_rna, bulge_type
   *
   * CFD score is intentionally absent; the Python scoring layer appends
   * it as the final column after computing it over the returned rows.
   *
   * @param sep  Field separator character (default: @c '\\t').
   * @return     A single row string with no trailing newline.
   */
  [[nodiscard]] std::string to_tsv_row(char sep = '\t') const;

  /**
   * @brief Column header names matching the @c to_tsv_row() column order.
   *
   * @return Vector of column name strings.
   */
  [[nodiscard]] static std::vector<std::string> tsv_header();

  // =====================================================================
  // Debug / diagnostic
  // =====================================================================

  /**
   * @brief Compact locus string for log messages.
   *
   * @return E.g. @c "chr1:123456(+)".
   */
  [[nodiscard]] std::string locus_string() const;

private:
  std::string chrom_;
  int32_t pos_;
  Strand strand_;
  std::string grna_;
  std::string target_;
  int mismatches_;
  int bulge_dna_;
  int bulge_rna_;
};

} // namespace crispritz

// =============================================================================
// std::hash specialisation (enables std::unordered_set<OffTarget>)
// =============================================================================

namespace std {
template <> struct hash<crispritz::OffTarget> {
  /**
   * @brief Hash over the locus identity fields that define equality.
   *
   * Uses the Boost hash_combine pattern to mix individual field hashes.
   * Hash collisions are expected (pigeonhole principle) but rare with
   * this mixing function for typical genomic data distributions.
   */
  std::size_t operator()(const crispritz::OffTarget &ot) const noexcept {
    // Boost-style hash_combine lambda
    auto combine = [](std::size_t seed, std::size_t v) noexcept -> std::size_t {
      return seed ^ (v + 0x9e3779b9u + (seed << 6) + (seed >> 2));
    };

    std::size_t h = std::hash<std::string>{}(ot.chrom());
    h = combine(h, std::hash<int32_t>{}(ot.pos()));
    h = combine(h, std::hash<char>{}(crispritz::to_char(ot.strand())));
    h = combine(h, std::hash<std::string>{}(ot.target()));
    return h;
  }
};
} // namespace std