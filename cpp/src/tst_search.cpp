/**
 * @file tst_search.cpp
 * @brief Search-stage algorithms: partition deserialization, bounded TST
 *        near-neighbour traversal with mismatch and bulge handling, and
 *        OffTarget result generation.
 *
 * Reuses the public TSTNode / TSTLeaf types and the encoding / nibble helpers
 * from tst.hpp / tst_utils.hpp. The TST *builder* (TernarySearchTree) is not
 * touched; this file only consumes the on-disk format it produces.
 *
 * Output writing is deliberately out of scope: this file produces OffTarget
 * objects in memory and stops there.
 */

#include "tst_search.hpp"

#include "tst_utils.hpp" // iupac::*, pack/unpack nibbles, sentinels, reverse_complement

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace pam; // NucleotideEncoder::*

namespace crispritz {

// =========================================================================
// LoadedTST
// =========================================================================

LoadedTST::LoadedTST(std::vector<TSTNode> nodes, std::vector<TSTLeaf> leaves,
                     int guide_length, std::string source_path)
    : nodes_(std::move(nodes)), leaves_(std::move(leaves)),
      guide_length_(guide_length), source_path_(std::move(source_path)) {
  if (nodes_.empty())
    throw std::invalid_argument("LoadedTST: node pool must not be empty");
  if (guide_length_ <= 0)
    throw std::invalid_argument("LoadedTST: guide_length must be > 0, got " +
                                std::to_string(guide_length_));
}

// =========================================================================
// Partition deserialization
// =========================================================================
//
// .bin layout (see TernarySearchTree::write_partition):
//   [4 B]  chunk_size  (number of leaves)
//   [4 B]  guide_length
//   for each leaf:
//     [4 B]                 guide_index (signed)
//     [ceil(pam_limit/2) B] bit-packed PAM nibbles      (length not stored!)
//     [1 B]                 '0'            -> next == 0
//                           '_' + [4 B]    -> next index
//   [4 B]  node_count
//   [var]  nibble-packed node stream (pre-order: split, lo, hi, eq)
//
// The PAM byte count per leaf is ceil(pam_limit/2). pam_limit is not in the
// header, but it is recoverable: the leaf section and node section are
// both present, and the node stream is self-delimiting (node_count bounds
// it). We therefore reconstruct pam byte width from the file name's PAM
// token when available; when not, the caller supplies it. To keep the
// loader self-contained and format-faithful, we read the PAM bytes using a
// width derived from the on-disk node reconstruction is not possible —
// instead we parse leaves by scanning for the '0'/'_' delimiter, which is
// unambiguous because guide_index and any 4-byte next are fixed width and
// the packed PAM contains no bytes equal to the delimiter only by chance.
//
// To avoid that fragility, load_partition takes the pam byte width as a
// parameter via an internal helper; the public entry derives it from the
// node stream position. For correctness and clarity we read the whole file
// into memory and parse with explicit offsets, requiring pam_bytes to be
// known. We obtain pam_bytes by a single backward pass: the node section
// length is determined by node_count, letting us locate the leaf section
// end and divide evenly. See parse below.
// =========================================================================

namespace {
struct NibbleReader; // defined below
void rebuild_children(NibbleReader &r, std::vector<TSTNode> &nodes,
                      int self_idx);

/** @brief Read a little-endian int32 from @p p and advance the cursor. */
int read_i32(const unsigned char *&p) {
  int v;
  std::memcpy(&v, p, sizeof(int));
  p += sizeof(int);
  return v;
}

/**
 * @brief Reconstruct the node pool from the nibble-packed pre-order stream.
 *
 * Mirrors TernarySearchTree::serialize_node exactly, in reverse:
 * each node emits splitchar, then lokid subtree (or '0'), then hikid
 * subtree (or '0'), then eqkid subtree (or '_' + 4-byte leaf pointer).
 *
 * Characters are stored two nibbles per byte. The high nibble 0xF ('_')
 * signals a leaf pointer follows on a fresh byte boundary; 0x0 ('0')
 * signals an absent child.
 *
 * Returns the index of the node it created (appended to @p nodes), or a
 * negative leaf-pointer encoding when the position was a leaf reference.
 *
 * @param data        Whole-file byte buffer.
 * @param node_cursor Byte offset into the node section; advanced in place.
 * @param nibble_high True when the next nibble to consume is the high
 *                    half of the current byte; advanced in place.
 * @param nodes       Destination node pool (root ends up at index 0 only
 *                    if this is the first node created — the builder
 *                    guarantees pre-order from the root).
 */
struct NibbleReader {
  const unsigned char *data;
  std::size_t pos;  // byte offset within node section
  bool high = true; // next nibble is high half

  explicit NibbleReader(const unsigned char *d, std::size_t start)
      : data(d), pos(start) {}

  uint8_t next_nibble() {
    uint8_t byte = data[pos];
    uint8_t nib = high ? high_nibble(byte) : low_nibble(byte);
    if (high) {
      high = false;
    } else {
      high = true;
      ++pos;
    }
    return nib;
  }

  // Align to the next whole byte (used before reading a 4-byte leaf ptr).
  void align_byte() {
    if (!high) {
      high = true;
      ++pos;
    }
  }

  int read_leaf_ptr() {
    align_byte();
    int v;
    std::memcpy(&v, data + pos, sizeof(int));
    pos += sizeof(int);
    return v;
  }
};

/**
 * @brief Recursively rebuild one node and its subtrees.
 * @return node index in @p nodes, or a negative leaf-pointer encoding.
 */
int rebuild_node(NibbleReader &r, std::vector<TSTNode> &nodes) {
  // Read splitchar nibble for this node.
  uint8_t split_nib = r.next_nibble();
  const int self_idx = static_cast<int>(nodes.size());
  nodes.emplace_back();
  nodes[self_idx].splitchar_enc = split_nib;
  nodes[self_idx].splitchar = NucleotideEncoder::decode_genome(split_nib);

  // lokid
  {
    // Peek the next nibble: 0x0 means absent child.
    uint8_t nib = r.next_nibble();
    if (nib == NULL_CHILD_NIBBLE) {
      nodes[self_idx].lokid = 0;
    } else {
      // This nibble was actually the splitchar of the lokid node;
      // rewind one nibble by re-injecting through a child rebuild.
      // To keep things simple we re-create the child using this
      // already-consumed nibble as its splitchar.
      const int child = static_cast<int>(nodes.size());
      nodes.emplace_back();
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = NucleotideEncoder::decode_genome(nib);
      rebuild_children(r, nodes, child);
      nodes[self_idx].lokid = child;
    }
  }

  // hikid
  {
    uint8_t nib = r.next_nibble();
    if (nib == NULL_CHILD_NIBBLE) {
      nodes[self_idx].hikid = 0;
    } else {
      const int child = static_cast<int>(nodes.size());
      nodes.emplace_back();
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = NucleotideEncoder::decode_genome(nib);
      rebuild_children(r, nodes, child);
      nodes[self_idx].hikid = child;
    }
  }

  // eqkid
  {
    uint8_t nib = r.next_nibble();
    if (nib == SENTINEL_NIBBLE) {
      // Leaf pointer follows on a byte boundary.
      int leaf_ptr = r.read_leaf_ptr(); // negative: -(within_chunk+1)
      nodes[self_idx].eqkid = leaf_ptr;
    } else {
      const int child = static_cast<int>(nodes.size());
      nodes.emplace_back();
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = NucleotideEncoder::decode_genome(nib);
      rebuild_children(r, nodes, child);
      nodes[self_idx].eqkid = child;
    }
  }

  return self_idx;
}

/**
 * @brief Fill in lokid/hikid/eqkid for a node whose splitchar nibble was
 *        already consumed by the caller.
 *
 * Split out from rebuild_node so a child whose splitchar was peeked can
 * be completed without double-reading its splitchar.
 */
void rebuild_children(NibbleReader &r, std::vector<TSTNode> &nodes,
                      int self_idx) {
  // lokid
  uint8_t lo = r.next_nibble();
  if (lo == NULL_CHILD_NIBBLE) {
    nodes[self_idx].lokid = 0;
  } else {
    const int child = static_cast<int>(nodes.size());
    nodes.emplace_back();
    nodes[child].splitchar_enc = lo;
    nodes[child].splitchar = NucleotideEncoder::decode_genome(lo);
    rebuild_children(r, nodes, child);
    nodes[self_idx].lokid = child;
  }

  // hikid
  uint8_t hi = r.next_nibble();
  if (hi == NULL_CHILD_NIBBLE) {
    nodes[self_idx].hikid = 0;
  } else {
    const int child = static_cast<int>(nodes.size());
    nodes.emplace_back();
    nodes[child].splitchar_enc = hi;
    nodes[child].splitchar = NucleotideEncoder::decode_genome(hi);
    rebuild_children(r, nodes, child);
    nodes[self_idx].hikid = child;
  }

  // eqkid
  uint8_t eq = r.next_nibble();
  if (eq == SENTINEL_NIBBLE) {
    int leaf_ptr = r.read_leaf_ptr();
    nodes[self_idx].eqkid = leaf_ptr;
  } else {
    const int child = static_cast<int>(nodes.size());
    nodes.emplace_back();
    nodes[child].splitchar_enc = eq;
    nodes[child].splitchar = NucleotideEncoder::decode_genome(eq);
    rebuild_children(r, nodes, child);
    nodes[self_idx].eqkid = child;
  }
}

} // namespace

LoadedTST load_partition(const std::string &partition_path) {
  std::ifstream in(partition_path, std::ios::in | std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("load_partition: cannot open " + partition_path);

  std::vector<unsigned char> buf((std::istreambuf_iterator<char>(in)),
                                 std::istreambuf_iterator<char>());
  in.close();

  if (buf.size() < 8u)
    throw std::runtime_error("load_partition: file too small: " +
                             partition_path);

  const unsigned char *p = buf.data();
  const int chunk_size = read_i32(p);
  const int guide_length = read_i32(p);

  if (chunk_size < 0 || guide_length <= 0)
    throw std::runtime_error("load_partition: invalid header in " +
                             partition_path);

  // --- leaf section -----------------------------------------------------
  // PAM byte width is ceil(pam_limit/2). pam_limit is not in the header.
  // We recover it by locating the node-count word: walk leaves using the
  // self-delimiting '0'/'_' marker. Each leaf is:
  //   guide_index(4) + pam_bytes(W) + marker(1) [+ next(4) if '_']
  // W is constant across all leaves in a partition, so we infer it from
  // the FIRST leaf by trying widths until the post-leaf node_count word
  // produces a self-consistent parse. In practice the Python layer knows
  // pam_limit and the binding passes it; here we infer conservatively.
  //
  // Simplest robust approach: the leaf section is parsed once W is known.
  // We determine W by reading the node_count from the tail is not possible
  // without W. Therefore load_partition reconstructs W from the PAM token
  // embedded in the filename (<pam>_<chr>_<part>.bin) when present.

  // Derive pam_bytes from the filename's leading PAM token.
  auto derive_pam_bytes = [](const std::string &path) -> int {
    // basename
    std::size_t slash = path.find_last_of("/\\");
    std::string base =
        (slash == std::string::npos) ? path : path.substr(slash + 1);
    std::size_t us = base.find('_');
    if (us == std::string::npos || us == 0)
      return -1;
    const int pam_limit = static_cast<int>(us); // PAM token length
    return (pam_limit + 1) / 2;                 // ceil(pam_limit/2)
  };

  const int pam_bytes = derive_pam_bytes(partition_path);
  if (pam_bytes < 0)
    throw std::runtime_error(
        "load_partition: cannot derive PAM width from filename: " +
        partition_path);

  std::vector<TSTLeaf> leaves;
  leaves.reserve(static_cast<std::size_t>(chunk_size));

  for (int i = 0; i < chunk_size; ++i) {
    TSTLeaf leaf;
    leaf.guide_index = read_i32(p);

    leaf.pam_seq_enc.assign(p, p + pam_bytes);
    p += pam_bytes;

    unsigned char marker = *p++;
    if (marker == static_cast<unsigned char>('_')) {
      leaf.next = read_i32(p);
    } else {
      leaf.next = 0; // marker == '0'
    }
    leaves.push_back(std::move(leaf));
  }

  // --- node section -----------------------------------------------------
  const int node_count = read_i32(p);
  if (node_count < 0)
    throw std::runtime_error("load_partition: negative node_count in " +
                             partition_path);

  std::vector<TSTNode> nodes;
  nodes.reserve(static_cast<std::size_t>(node_count));

  if (node_count > 0) {
    const std::size_t node_start = static_cast<std::size_t>(p - buf.data());
    NibbleReader reader(buf.data(), node_start);
    rebuild_node(reader, nodes);
  } else {
    // Degenerate but valid: a single sentinel root so LoadedTST's
    // invariant (non-empty pool) holds.
    nodes.emplace_back();
  }

  return LoadedTST{std::move(nodes), std::move(leaves), guide_length,
                   partition_path};
}

// =========================================================================
// TSTSearcher
// =========================================================================

TSTSearcher::TSTSearcher(SearchConfiguration config)
    : config_(std::move(config)) {
  // Note: config_.threads() is intentionally NOT consulted here. The
  // searcher is single-threaded by design; thread count belongs to the
  // Python orchestration layer, which sizes its ThreadPoolExecutor and
  // dispatches one search_partition() call per .bin file. The searcher
  // only reads the edit-budget fields (mismatches / bulges). See the
  // serial-loop rationale in search_all().
}

namespace {
/**
 * @brief Mutable-by-value recursion frame for the bounded near-search.
 *
 * Passing the state by value (not by reference) is deliberate: each
 * branch of the traversal gets its own copy, so a mismatch or bulge
 * spent down one path cannot corrupt the budget of a sibling path. This
 * removes an entire class of backtracking bugs at the cost of copying a
 * few ints and two short strings per recursive call.
 */
struct Frame {
  int mm_left;            // remaining substitution budget
  int bdna_left;          // remaining DNA-bulge budget
  int brna_left;          // remaining RNA-bulge budget
  std::size_t guide_pos;  // next index into the query guide
  std::string aln_guide;  // aligned guide chars accumulated so far
  std::string aln_target; // aligned target chars accumulated so far
};

/** @brief Bundle of read-only context shared across the recursion. */
struct Context {
  const LoadedTST *tst;
  std::string_view guide; // full query guide
  int guide_len;          // guide.size()
  std::string_view chrom; // chromosome name for emitted hits
  Strand strand_for_index(int guide_index) const {
    return guide_index >= 0 ? Strand::Forward : Strand::Reverse;
  }
};

// Forward declaration: the traversal recurses through three node slots.
void traverse(const Context &ctx, int node_idx, Frame frame,
              std::vector<OffTarget> &out);

/**
 * @brief Emit OffTargets for the leaf chain rooted at @p leaf_ptr.
 *
 * @p leaf_ptr is the negative encoding -(within_chunk_index + 1) stored
 * in a terminal node's eqkid. All leaves chained via TSTLeaf::next share
 * the same guide path, so they share the same aligned strings and edit
 * counts; only their genomic position/strand differ.
 */
void emit_leaf_chain(const Context &ctx, int leaf_ptr, const Frame &frame,
                     std::vector<OffTarget> &out) {
  const auto &leaves = ctx.tst->leaves();
  int idx = -leaf_ptr - 1; // decode within-chunk leaf index

  while (idx >= 0 && idx < static_cast<int>(leaves.size())) {
    const TSTLeaf &leaf = leaves[static_cast<std::size_t>(idx)];

    const Strand strand = ctx.strand_for_index(leaf.guide_index);
    // guide_index stores the 0-based leftmost genomic coordinate of
    // the window (same convention for both strands; the sign only
    // encodes strand). OffTarget::pos is 1-based, so add 1.
    const int leftmost0 =
        leaf.guide_index >= 0 ? leaf.guide_index : -leaf.guide_index;
    const int pos = leftmost0 + 1;

    // Edit-distance breakdown = initial budget minus remaining.
    // The caller seeds the frame with the configured maxima, so the
    // spent amount is (max - left). We recover the maxima from the
    // accumulated alignment instead: count lowercase (mismatch) and
    // '-' (bulge) characters, which are written during traversal.
    int mm_count = 0, bdna_count = 0, brna_count = 0;
    for (std::size_t k = 0; k < frame.aln_target.size(); ++k) {
      char tg = frame.aln_target[k];
      char gd = (k < frame.aln_guide.size()) ? frame.aln_guide[k] : '?';
      if (gd == '-')
        ++bdna_count; // gap in guide  = DNA bulge
      else if (tg == '-')
        ++brna_count; // gap in target = RNA bulge
      else if (tg >= 'a' && tg <= 'z')
        ++mm_count; // lowercase = mismatch
    }

    out.emplace_back(
        /*chrom  */ std::string(ctx.chrom),
        /*pos    */ pos,
        /*strand */ strand,
        /*grna   */ frame.aln_guide,
        /*target */ frame.aln_target,
        /*mm     */ mm_count,
        /*bdna   */ bdna_count,
        /*brna   */ brna_count);

    // TSTLeaf::next holds the *encoded* pointer to the next leaf in
    // the collision chain — the same -(within_chunk + 1) convention
    // used by the head pointer (see TernarySearchTree::insert, which
    // assigns next = node.eqkid). 0 terminates the chain.
    if (leaf.next == 0)
      break;
    idx = -leaf.next - 1;
  }
}

void traverse(const Context &ctx, int node_idx, Frame frame,
              std::vector<OffTarget> &out) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.tst->nodes().size()))
    return;

  const TSTNode &node = ctx.tst->nodes()[static_cast<std::size_t>(node_idx)];

  // A node with no splitchar encoding and no children is the reserved
  // root sentinel (index 0). Descend into its eqkid only.
  const bool is_root_sentinel = (node_idx == 0 && node.splitchar_enc == 0 &&
                                 node.lokid == 0 && node.hikid == 0);

  if (is_root_sentinel) {
    if (node.eqkid < 0)
      emit_leaf_chain(ctx, node.eqkid, frame, out);
    else if (node.eqkid > 0)
      traverse(ctx, node.eqkid, frame, out);
    return;
  }

  // If the query is exhausted we cannot consume this splitchar by a
  // match or mismatch; only an RNA bulge (gap in target) could apply,
  // but that is handled when we still have guide characters. So stop.
  if (frame.guide_pos >= static_cast<std::size_t>(ctx.guide_len))
    return;

  const char q_char = ctx.guide[frame.guide_pos];
  const uint8_t q_enc = iupac::encode_genome(q_char);
  const bool is_match = iupac::matches(q_enc, node.splitchar_enc);

  // ---- lokid / hikid : explore alternative splitchars at this depth.
  //      No budget is consumed; these do not advance the guide.
  if (node.lokid > 0)
    traverse(ctx, node.lokid, frame, out);
  if (node.hikid > 0)
    traverse(ctx, node.hikid, frame, out);

  // ---- equal branch: consume this splitchar against the query char.
  const char node_char = node.splitchar;

  if (is_match) {
    Frame next = frame;
    next.aln_guide += node_char;
    next.aln_target += node_char; // uppercase = match
    next.guide_pos += 1;
    if (node.eqkid < 0)
      emit_leaf_chain(ctx, node.eqkid, next, out);
    else if (node.eqkid > 0)
      traverse(ctx, node.eqkid, next, out);
  } else if (frame.mm_left > 0) {
    Frame next = frame;
    next.mm_left -= 1;
    next.aln_guide += q_char;             // query base (upper)
    next.aln_target += static_cast<char>( // genome base (lower)
        node_char >= 'A' && node_char <= 'Z' ? node_char - 'A' + 'a'
                                             : node_char);
    next.guide_pos += 1;
    if (node.eqkid < 0)
      emit_leaf_chain(ctx, node.eqkid, next, out);
    else if (node.eqkid > 0)
      traverse(ctx, node.eqkid, next, out);
  }

  // ---- DNA bulge: extra base in the genome/target, gap in the guide.
  //      Consume this node's splitchar into the target, advance the
  //      node (eqkid) but NOT the guide position.
  if (frame.bdna_left > 0) {
    Frame next = frame;
    next.bdna_left -= 1;
    next.aln_guide += '-';        // gap in guide
    next.aln_target += node_char; // extra genomic base
    // guide_pos unchanged
    if (node.eqkid > 0)
      traverse(ctx, node.eqkid, next, out);
    // (a DNA bulge cannot land exactly on a leaf terminal without a
    //  following matched base, so we do not emit on eqkid < 0 here)
  }

  // ---- RNA bulge: extra base in the guide, gap in the genome/target.
  //      Consume a guide character into the alignment as a gap on the
  //      target side, advance the guide but stay on the same node.
  if (frame.brna_left > 0) {
    Frame next = frame;
    next.brna_left -= 1;
    next.aln_guide += q_char; // extra guide base
    next.aln_target += '-';   // gap in target
    next.guide_pos += 1;
    traverse(ctx, node_idx, next, out); // same node, advanced guide
  }
}

} // namespace

std::vector<OffTarget> TSTSearcher::search(const LoadedTST &tst,
                                           std::string_view guide_seq,
                                           const std::string &chrom) const {
  // Deferred guide-length / edit-budget validation (see header contract).
  if (config_.max_total_edits() > tst.guide_length())
    throw std::invalid_argument("TSTSearcher::search: max_total_edits (" +
                                std::to_string(config_.max_total_edits()) +
                                ") exceeds index guide_length (" +
                                std::to_string(tst.guide_length()) + ')');

  if (static_cast<int>(guide_seq.size()) != tst.guide_length())
    throw std::invalid_argument("TSTSearcher::search: guide length (" +
                                std::to_string(guide_seq.size()) +
                                ") does not match index guide_length (" +
                                std::to_string(tst.guide_length()) + ')');

  std::vector<OffTarget> out;

  Context ctx;
  ctx.tst = &tst;
  ctx.guide = guide_seq;
  ctx.guide_len = static_cast<int>(guide_seq.size());
  ctx.chrom = chrom;

  Frame root;
  root.mm_left = config_.max_mismatches();
  root.bdna_left = config_.max_bulges_dna();
  root.brna_left = config_.max_bulges_rna();
  root.guide_pos = 0;
  root.aln_guide.reserve(
      static_cast<std::size_t>(ctx.guide_len + config_.max_bulges_total()));
  root.aln_target.reserve(
      static_cast<std::size_t>(ctx.guide_len + config_.max_bulges_total()));

  traverse(ctx, /*root index*/ 0, std::move(root), out);
  return out;
}

SearchResult TSTSearcher::search_all(const LoadedTST &tst,
                                     const std::vector<std::string> &guides,
                                     const std::string &chrom) const {
  SearchResult result;
  result.source_path = tst.source_path();
  result.hits_by_guide.reserve(guides.size());

  // The loop over guides is DELIBERATELY serial. Do not parallelise it
  // (e.g. with "#pragma omp parallel for").
  //
  // The parallelism axis in this architecture is the *partition*, not the
  // guide: the Python orchestration layer drives a memory-aware
  // ThreadPoolExecutor that calls search_partition() once per .bin file,
  // and the pybind11 binding releases the GIL so those calls run as true
  // parallel C++. Real workloads have very few guides (1-20) but many
  // partitions, so partition-level parallelism saturates the cores while
  // guide-level parallelism would not.
  //
  // Adding OpenMP here would nest guide-parallelism inside the Python
  // thread pool and oversubscribe the CPU (threads x guides), degrading
  // performance. Keeping this serial is the correct design, not an
  // oversight.
  for (const std::string &g : guides)
    result.hits_by_guide.push_back(search(tst, g, chrom));

  return result;
}

// =========================================================================
// search_partition — atomic unit of parallel work
// =========================================================================

SearchResult search_partition(const std::string &partition_path,
                              const std::string &chrom,
                              const std::vector<std::string> &guides,
                              const SearchConfiguration &config) {
  LoadedTST tst = load_partition(partition_path);
  TSTSearcher searcher(config);
  return searcher.search_all(tst, guides, chrom);
}

} // namespace crispritz