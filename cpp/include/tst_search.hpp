#pragma once

#include "offtarget.hpp"            // crispritz::OffTarget
#include "search_configuration.hpp" // crispritz::SearchConfiguration
#include "tst.hpp"                  // crispritz::TSTNode, crispritz::TSTLeaf

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace crispritz {

// =========================================================================
// LoadedTST
// =========================================================================

/**
 * @brief An immutable, in-memory Ternary Search Tree loaded from one
 *        serialized @c .bin partition.
 *
 * The builder class @c TernarySearchTree (tst.hpp) only *writes* partitions;
 * it keeps its node/leaf arrays private and exposes no query path. The search
 * stage therefore needs its own read-only representation, reconstructed from
 * the same public @c TSTNode and @c TSTLeaf types so the two stages share one
 * definition of the on-disk layout.
 *
 * A @c LoadedTST owns its node and leaf arrays and is immutable after
 * construction: every accessor is @c const and there are no mutators. This is
 * what makes a single instance safe to read concurrently from many threads,
 * which is the property the partition-parallel search relies on (Python owns
 * the outer thread pool; each worker reads one @c LoadedTST).
 *
 * ## Index geometry travels with the index
 * @c guide_length() is read from the partition header at load time — it is an
 * indexing-time property, not a search parameter (see SearchConfiguration's
 * documentation). The search layer obtains it here rather than from the user.
 *
 * @note This class declares interface only; construction is performed by
 *       @c load_partition() (declared below). The constructor is intentionally
 *       not specified here beyond what the loader requires.
 */
class LoadedTST {
public:
  /**
   * @brief Construct from already-deserialized arrays and header values.
   *
   * Ownership of both vectors is transferred into the object.
   *
   * @param nodes         TST node pool; index 0 is the root.
   * @param leaves        Leaf array indexed by the negative @c eqkid
   *                      encoding used in @c TSTNode.
   * @param guide_length  Guide length read from the partition header (> 0).
   * @param source_path   Originating @c .bin path, retained for diagnostics.
   *
   * @throws std::invalid_argument if @p nodes is empty or
   *         @p guide_length <= 0.
   */
  LoadedTST(std::vector<TSTNode> nodes, std::vector<TSTLeaf> leaves,
            int guide_length, std::string source_path);

  /** @return The TST node pool (root at index 0). */
  [[nodiscard]] const std::vector<TSTNode> &nodes() const noexcept {
    return nodes_;
  }

  /** @return The leaf array. */
  [[nodiscard]] const std::vector<TSTLeaf> &leaves() const noexcept {
    return leaves_;
  }

  /** @return Guide length recorded in the partition header (> 0). */
  [[nodiscard]] int guide_length() const noexcept { return guide_length_; }

  /** @return Originating @c .bin file path (for diagnostics/logging). */
  [[nodiscard]] const std::string &source_path() const noexcept {
    return source_path_;
  }

  /** @return Number of nodes in the pool. */
  [[nodiscard]] std::size_t node_count() const noexcept {
    return nodes_.size();
  }

  /** @return Number of leaves in the array. */
  [[nodiscard]] std::size_t leaf_count() const noexcept {
    return leaves_.size();
  }

private:
  std::vector<TSTNode> nodes_;
  std::vector<TSTLeaf> leaves_;
  int guide_length_;
  std::string source_path_;
};

// =========================================================================
// load_partition
// =========================================================================

/**
 * @brief Deserialize one @c .bin partition into a @c LoadedTST.
 *
 * Inverts the serialization performed by @c TernarySearchTree::save():
 * reads the header (leaf count, guide length), the leaf array, the node
 * count, and the nibble-packed node stream.
 *
 * @param partition_path  Path to a single @c .bin partition file.
 * @return                A populated, immutable @c LoadedTST.
 *
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 *
 * @note Interface only — no implementation is provided in this phase.
 */
[[nodiscard]] LoadedTST load_partition(const std::string &partition_path);

// =========================================================================
// SearchResult
// =========================================================================

/**
 * @brief Off-target hits for all guides against a single partition, plus
 *        lightweight diagnostics.
 *
 * Returning a bare @c std::vector<OffTarget> would lose three things the
 * caller needs in the multi-guide, multi-partition setting:
 *   - which guide each hit belongs to (hits are grouped per guide),
 *   - the originating partition (for logging and provenance), and
 *   - a cheap total-hit count without re-summing the nested vectors.
 *
 * @c hits_by_guide is indexed in lockstep with the @c guides argument passed
 * to the search: @c hits_by_guide[i] holds every @c OffTarget found for
 * @c guides[i] in this partition (possibly empty).
 */
struct SearchResult {
  /// Per-guide hit lists; outer index matches the input guide order.
  std::vector<std::vector<OffTarget>> hits_by_guide;

  /// Path of the partition these hits came from (provenance/logging).
  std::string source_path;

  /**
   * @brief Total number of hits across all guides in this partition.
   * @return Sum of @c hits_by_guide[i].size() over all i.
   */
  [[nodiscard]] std::size_t total_hits() const noexcept {
    std::size_t n = 0;
    for (const auto &g : hits_by_guide)
      n += g.size();
    return n;
  }

  /**
   * @brief Number of guides this result covers.
   * @return hits_by_guide.size().
   */
  [[nodiscard]] std::size_t guide_count() const noexcept {
    return hits_by_guide.size();
  }
};

// =========================================================================
// TSTSearcher
// =========================================================================

/**
 * @brief Stateless near-neighbour searcher over an already-loaded TST.
 *
 * Separated from the partition entry point so the traversal can be unit
 * tested against an in-memory @c LoadedTST without any filesystem access.
 * The searcher holds only the @c SearchConfiguration it was constructed with;
 * it allocates no persistent state, so a single instance is reentrant and may
 * be used by multiple threads, each searching a different @c LoadedTST.
 *
 * ## Edit-budget vs guide-length check
 * The configuration alone cannot validate that @c max_total_edits() fits
 * within the guide length, because the guide length lives in the index. That
 * check is performed here, against @c LoadedTST::guide_length(), before the
 * traversal begins (see @c search()).
 */
class TSTSearcher {
public:
  /**
   * @brief Construct a searcher bound to a configuration.
   * @param config  Validated search parameters (copied in).
   */
  explicit TSTSearcher(SearchConfiguration config);

  /**
   * @brief Search one loaded partition for one guide.
   *
   * @param tst        The loaded, read-only index to search.
   * @param guide_seq  The query guide in the same canonical orientation
   *                   used during indexing (see @c TSTLeaf::guide_seq).
   * @return           Every off-target for this guide in this partition.
   *
   * @throws std::invalid_argument if @c config.max_total_edits() exceeds
   *         @c tst.guide_length(), or if @p guide_seq length is
   *         inconsistent with @c tst.guide_length().
   *
   * @param tst        The loaded, read-only index to search.
   * @param guide_seq  Query guide in canonical orientation.
   * @param chrom      Chromosome name to record in each emitted OffTarget.
   *                   Resolved by the Python layer from the partition
   *                   filename and passed down, so results need no
   *                   post-hoc chromosome fixup.
   */
  [[nodiscard]] std::vector<OffTarget> search(const LoadedTST &tst,
                                              std::string_view guide_seq,
                                              const std::string &chrom) const;

  /**
   * @brief Search one loaded partition for many guides.
   *
   * Convenience wrapper that applies @c search() to each guide and groups
   * the results. The returned @c SearchResult is indexed in lockstep with
   * @p guides.
   *
   * @param tst     The loaded, read-only index to search.
   * @param guides  Query guides in canonical orientation.
   * @param chrom   Chromosome name recorded in every emitted OffTarget.
   * @return        Per-guide hit lists plus partition provenance.
   *
   * @throws std::invalid_argument under the same conditions as @c search().
   */
  [[nodiscard]] SearchResult search_all(const LoadedTST &tst,
                                        const std::vector<std::string> &guides,
                                        const std::string &chrom) const;

  /** @return The configuration this searcher was constructed with. */
  [[nodiscard]] const SearchConfiguration &config() const noexcept {
    return config_;
  }

private:
  SearchConfiguration config_;
};

// =========================================================================
// search_partition — atomic unit of parallel work
// =========================================================================

/**
 * @brief Load one partition and search it for all guides, in one call.
 *
 * This free function is the single primitive the Python orchestration layer
 * invokes per partition. Python owns the outer concurrency (a memory-aware
 * @c ThreadPoolExecutor), calling @c search_partition() once per @c .bin file
 * across worker threads; the pybind11 binding releases the GIL so these calls
 * run as true parallel C++. Each call is fully independent — it loads its own
 * @c LoadedTST and shares no mutable state with any other call — so no locking
 * is required.
 *
 * Equivalent to:
 * @code
 *   LoadedTST   tst = load_partition(partition_path);
 *   TSTSearcher searcher(config);
 *   return searcher.search_all(tst, guides);
 * @endcode
 * but exposed as one entry point so the binding surface stays minimal.
 *
 * @param partition_path  Path to a single @c .bin partition file.
 * @param chrom           Chromosome name (resolved Python-side from the
 *                        partition filename) recorded in every OffTarget.
 * @param guides          Query guides in canonical orientation.
 * @param config          Validated search parameters.
 * @return                Per-guide hit lists plus partition provenance.
 *
 * @throws std::runtime_error    if the partition cannot be loaded.
 * @throws std::invalid_argument if the edit budget is incompatible with the
 *                               loaded guide length.
 */
[[nodiscard]] SearchResult
search_partition(const std::string &partition_path, const std::string &chrom,
                 const std::vector<std::string> &guides,
                 const SearchConfiguration &config);

} // namespace crispritz