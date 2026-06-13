#pragma once

#include "tst_utils.hpp"

#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

namespace crispritz {

// =========================================================================
// TSTNode
// =========================================================================

/**
 * @brief A single internal node of the Ternary Search Tree.
 *
 * Indices are used instead of raw pointers so that the node array can be
 * stored in a contiguous std::vector without invalidating references on
 * reallocation.  A value of 0 in any child field means "no child" (the root
 * is always at index 0, so child 0 is reserved as the null sentinel).
 *
 * The eqkid field doubles as a leaf pointer: when negative, the node IS a
 * leaf and the leaf index is (-eqkid - 1), matching the legacy encoding
 * used by the binary .bin format and searchOnTST.cpp.
 *
 * @note splitchar_enc stores the 4-bit IUPAC encoding (iupac::encode_genome)
 *       of the node character so that match checks during search are a single
 *       bitwise AND, avoiding a second lookup.
 */
struct TSTNode {
  char splitchar = '\0';     ///< IUPAC character stored at this node.
  uint8_t splitchar_enc = 0; ///< 4-bit encoding of splitchar.
  int lokid = 0;             ///< Index of the left  (less-than) child.
  int hikid = 0;             ///< Index of the right (greater-than) child.
  int eqkid = 0;             ///< Index of the equal child, or -(leaf+1).
};

// =========================================================================
// TSTLeaf
// =========================================================================

/**
 * @brief Metadata stored at a TST leaf, one per PAM site extracted from the
 *        genome.
 *
 * The legacy code stored leaves in a parallel array of Tleaf structs indexed
 * via the negative eqkid values of terminal TSTNodes.  That layout is
 * preserved here so the binary serialization format does not change.
 *
 * @param guide_index  Genomic position of the target site.  Positive values
 *                     denote the forward strand; negative values denote the
 *                     reverse strand (following the convention in the legacy
 *                     code and searchOnTST.cpp).
 * @param guide_seq    The guide-length subsequence extracted from the genome,
 *                     stored 5'→3' AFTER reversal/complementation so that
 *                     lexicographic sort and TST insertion operate on the
 *                     canonical form used during search.
 * @param pam_seq_enc  Bit-packed PAM nucleotides (two per byte, high nibble
 *                     first) exactly as written to the .bin leaf section.
 * @param next         Index of the next leaf sharing the same guide_seq but a
 *                     different PAM site (linked-list of collisions). 0 means
 *                     no next leaf.
 */
struct TSTLeaf {
  int guide_index = 0;
  std::string guide_seq;            ///< Guide subsequence (canonical order).
  std::vector<uint8_t> pam_seq_enc; ///< Bit-packed PAM bytes.
  int next = 0;                     ///< Next leaf with same guide_seq.
};

// =========================================================================
// TernarySearchTree
// =========================================================================

/**
 * @brief Builds and serializes a Ternary Search Tree index from a genomic
 *        sequence and a PAM specification.
 *
 * Typical usage:
 * @code
 *   TernarySearchTree tst(sequence, pam, pam_length, pam_limit,
 *                         pam_at_start, outdir, max_bulges, num_threads);
 *   tst.build();   // PAM search + extraction + sort + insert + save
 * @endcode
 *
 * The serialized .bin files are consumed by the existing searchOnTST binary
 * without modification; the format is byte-for-byte identical to the legacy
 * mainParallel.cpp output.
 *
 * Thread safety: build() uses OpenMP to parallelize PAM search.
 * All other methods are single-threaded and must not be called concurrently.
 */
class TernarySearchTree {
public:
  /**
   * @brief Construct the TST builder.
   *
   * No computation is performed here; call build() to run the pipeline.
   *
   * @param sequence    Full genomic sequence (single chromosome, uppercase).
   * @param chr_name    Chromosome / contig identifier used in output filenames.
   * @param pam_seq     The PAM-only string (e.g. "NGG"), without guide Ns.
   * @param pam_length  Total length of the PAM+guide pattern (e.g. 23 for
   *                    NNNNNNNNNNNNNNNNNNNNNGG).
   * @param pam_limit   Length of the PAM portion only (e.g. 3 for "NGG").
   * @param pam_at_start  True when the PAM precedes the guide (e.g. Cas12a).
   * @param outdir      Directory where .bin partition files will be written.
   *                    The directory must already exist.
   * @param max_bulges  Maximum number of bulges; extra bases extracted per
   *                    site to allow bulge-aware search later.
   * @param num_threads Number of OpenMP threads for PAM search.
   *
   * @throws std::runtime_error if guide_length ≤ 0 or pam_limit ≤ 0 or
   *         outdir is empty.
   */
  TernarySearchTree(std::string_view sequence, std::string_view chr_name,
                    std::string_view pam_seq, int pam_length, int pam_limit,
                    bool pam_at_start, std::string_view outdir,
                    int max_bulges = 0, int num_threads = 1);

  /**
   * @brief Run the full build pipeline:
   *        1. PAM search  → site positions
   *        2. Sequence extraction → TSTLeaf array
   *        3. Lexicographic sort
   *        4. Balanced median insertion into TST nodes
   *        5. Serialization to .bin file(s) under outdir_
   *
   * @throws std::runtime_error if no PAM sites are found after N-filtering,
   *         or if a .bin file cannot be opened for writing.
   */
  void build();

  /**
   * @brief Write the constructed TST to disk.
   *
   * Large chromosomes are split into chunks of LEAVES_PER_GROUP leaves,
   * each written to a separate file named:
   * @code
   *   <outdir_>/<pam_seq_>_<chr_name_>_<part>.bin
   * @endcode
   *
   * This method is called internally by build() and is exposed for testing.
   *
   * @throws std::runtime_error on I/O failure.
   */
  void save() const;

  /** @return Number of valid leaves after extraction (excluding discarded Ns).
   */
  int leaf_count() const noexcept { return static_cast<int>(leaves_.size()); }

  /** @return Number of TST nodes allocated during the last build(). */
  int node_count() const noexcept { return static_cast<int>(nodes_.size()); }

  /** @brief Maximum leaves stored in a single .bin partition. */
  static constexpr int LEAVES_PER_GROUP = 5'000'000;

private:
  // ---------------------------------------------------------------------
  // Configuration (set in constructor, immutable after construction)
  // ---------------------------------------------------------------------
  std::string sequence_;
  std::string chr_name_;
  std::string pam_seq_;
  int pam_length_;
  int pam_limit_;
  int guide_length_; ///< pam_length_ - pam_limit_
  bool pam_at_start_;
  std::string outdir_; ///< Output directory for .bin files.
  int max_bulges_;
  int num_threads_;

  // ------------------------------------------------------------------
  // Build state (mutated by build() / save())
  // ------------------------------------------------------------------
  std::vector<TSTLeaf> leaves_; ///< Extracted and sorted target leaves.
  std::vector<TSTNode> nodes_;  ///< TST node pool (index 0 = root).
  int nodes_used_ = 0;

  // ------------------------------------------------------------------
  // Sequence extraction helpers
  // ------------------------------------------------------------------

  /**
   * @brief Extract a forward-strand target at genomic position @p pos.
   *
   * Reads pam_length_ + max_bulges_ characters from sequence_ at @p pos,
   * discards the window if any 'N' is present, splits into guide and PAM
   * subsequences, and appends a TSTLeaf to @p dest.
   *
   * @param pos   0-based genomic start position.
   * @param dest  Destination leaf vector.
   */
  void extract_forward(int pos, std::vector<TSTLeaf> &dest) const;

  /**
   * @brief Extract a reverse-strand target at genomic position @p pos.
   *
   * Reads pam_length_ + max_bulges_ characters, reverse-complements the
   * window, then applies the same split-and-push logic as extract_forward().
   *
   * @param pos   0-based genomic position (absolute value of the negative
   *              index returned by the PAM search).
   * @param dest  Destination leaf vector.
   */
  void extract_reverse(int pos, std::vector<TSTLeaf> &dest) const;

  /**
   * @brief Populate leaves_ from the PAM search result vector.
   *
   * Dispatches each hit to extract_forward() or extract_reverse() based on
   * its sign and the pam_at_start_ flag, matching the strand convention
   * used by mainParallel.cpp.
   *
   * @param pam_sites  Signed position vector from pam::search_pam_sites_fast.
   */
  void extract_sequences(const std::vector<int> &pam_sites);

  /**
   * @brief Encode a PAM substring into the bit-packed byte format used in
   *        TSTLeaf::pam_seq_enc and the .bin leaf section.
   *
   * Two IUPAC 4-bit codes are packed per byte (high nibble first).
   * An odd-length PAM occupies only the high nibble of the last byte.
   *
   * @param pam_str  Already-oriented PAM substring.
   * @return         Packed bytes of length ceil(pam_str.size() / 2).
   */
  static std::vector<uint8_t> encode_pam_bytes(std::string_view pam_str);

  // ---------------------------------------------------------------------
  // TST insertion helpers
  // ---------------------------------------------------------------------

  /**
   * @brief Allocate a fresh TST node and return its index.
   *
   * Appends a zero-initialised TSTNode to nodes_ and increments
   * nodes_used_.  The root (index 0) must be pre-allocated before the
   * first insert() call.
   *
   * @return Index of the newly allocated node.
   */
  int alloc_node();

  /**
   * @brief Insert the guide sequence of leaf @p leaf_idx into the TST.
   *
   * Navigates left/right on character mismatches and equal on match,
   * following the standard TST insertion algorithm.  On reaching the end
   * of the string the within-chunk leaf index is stored in eqkid using the
   * legacy encoding -(within_chunk_idx + 1).
   *
   * @param guide_str    Guide sequence to insert.
   * @param leaf_idx     Index into leaves_ of the leaf being inserted.
   * @param chunk_offset Index of the first leaf in the current partition
   *                     chunk (subtracted to get the within-chunk index).
   */
  void insert(std::string_view guide_str, int leaf_idx, int chunk_offset);

  /**
   * @brief Balanced median insertion of leaves[lo..hi] into the TST.
   *
   * Inserts the median first, then recurses on each half, producing a
   * balanced tree with O(log n) expected depth on sorted input.
   *
   * @param lo           First index of the range (inclusive).
   * @param hi           Last  index of the range (inclusive).
   * @param chunk_offset Start index of the current partition chunk.
   */
  void insert_balanced(int lo, int hi, int chunk_offset);

  // ---------------------------------------------------------------------
  // Serialization helpers
  // ---------------------------------------------------------------------

  /**
   * @brief Write one .bin partition to disk under outdir_.
   *
   * File layout (byte-exact match with legacy saveTST + serialize):
   * @code
   * [4 bytes]  number of leaves in this chunk
   * [4 bytes]  guide length (pam_length_ - pam_limit_)
   * for each leaf:
   *   [4 bytes]               guide_index (signed)
   *   [ceil(pam_limit_/2) B]  bit-packed PAM nucleotides
   *   [1 byte]                '0' if next==0, else '_' + [4 bytes] next
   * [4 bytes]  number of nodes
   * [variable] serialized TST nodes (writePair encoding)
   * @endcode
   *
   * @param part         1-based partition number (appended to filename).
   * @param chunk_start  Index of the first leaf in this partition.
   * @param chunk_end    One-past the last leaf in this partition.
   *
   * @throws std::runtime_error if the output file cannot be opened.
   */
  void write_partition(int part, int chunk_start, int chunk_end) const;

  /**
   * @brief Recursive pre-order TST serialization (mirrors legacy
   *        serialize()).
   *
   * Characters are buffered two-at-a-time (writePair equivalent) before
   * being flushed to the stream as a single packed byte.
   *
   * @param node_idx  Index of the node being serialized.
   * @param out       Binary output stream (must be open).
   * @param buf       Two-char accumulation buffer; pass {'\0','\0'} at root.
   * @param buf_pos   Current write position in buf (0 or 1); updated in place.
   */
  void serialize_node(int node_idx, std::ofstream &out, char buf[2],
                      int &buf_pos) const;

  /**
   * @brief Flush the two-char buffer as one packed byte and reset buf_pos.
   *
   * @param buf     Two-char buffer (buf[0] is high nibble, buf[1] is low).
   * @param buf_pos Must equal 2 on entry; reset to 0 on exit.
   * @param out     Binary output stream.
   */
  static void flush_pair(const char buf[2], int &buf_pos, std::ofstream &out);

  /**
   * @brief Buffer one character; flush when both slots are filled.
   *
   * @param c       Character to buffer ('0' for null child, '_' for sentinel).
   * @param buf     Two-char accumulation buffer.
   * @param buf_pos Current position (0 or 1); updated in place.
   * @param out     Binary output stream.
   */
  static void buffer_char(char c, char buf[2], int &buf_pos,
                          std::ofstream &out);
};

// =========================================================================
// Free function — pybind11 / Python API entry point
// =========================================================================

/**
 * @brief Build a TST index for a single genomic sequence and write the
 *        resulting .bin partition files to @p outdir.
 *
 * This is the function exposed through pybind11.  It constructs a
 * TernarySearchTree, calls build(), and lets save() write output under
 * @p outdir.
 *
 * @param sequence    Full genomic sequence string (single chromosome,
 * uppercase).
 * @param chr_name    Chromosome name (used in output filenames).
 * @param pam_seq     PAM-only string (e.g. "NGG").
 * @param pam_length  Full pattern length including guide Ns.
 * @param pam_limit   PAM-only length.
 * @param pam_at_start  True for PAM-upstream (e.g. Cas12a).
 * @param outdir      Directory where .bin files will be written (must exist).
 * @param max_bulges  Maximum bulge count (default 0).
 * @param num_threads Number of threads for PAM search (default 1).
 *
 * @throws std::runtime_error on any build or I/O failure.
 */
void build_tree(const std::string &sequence, const std::string &chr_name,
                const std::string &pam_seq, int pam_length, int pam_limit,
                bool pam_at_start, const std::string &outdir,
                int max_bulges = 0, int num_threads = 1);

} // namespace crispritz