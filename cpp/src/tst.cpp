#include "tst.hpp"
#include "pam_search.hpp"
#include "tst_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace crispritz {

// =========================================================================
// Constructor
// =========================================================================

TernarySearchTree::TernarySearchTree(std::string_view sequence,
                                     std::string_view chr_name,
                                     std::string_view pam_seq, int pam_length,
                                     int pam_limit, bool pam_at_start,
                                     std::string_view outdir, int max_bulges,
                                     int num_threads)
    : sequence_(sequence), chr_name_(chr_name), pam_seq_(pam_seq),
      pam_length_(pam_length), pam_limit_(pam_limit),
      guide_length_(pam_length - pam_limit), pam_at_start_(pam_at_start),
      outdir_(outdir), max_bulges_(max_bulges),
      num_threads_(num_threads > 0 ? num_threads : 1) {
  if (guide_length_ <= 0)
    throw std::runtime_error(
        "guide_length must be positive (pam_length > pam_limit required)");
  if (pam_limit_ <= 0)
    throw std::runtime_error("pam_limit must be positive");
  if (outdir_.empty())
    throw std::runtime_error("outdir must not be empty");
}

// =========================================================================
// Public entry point
// =========================================================================

void TernarySearchTree::build() {
  // --- 1. PAM search ---------------------------------------------------
  pam::SearchParams params(pam_length_, pam_limit_, pam_at_start_,
                           num_threads_);
  pam::CompactGenome genome_bits(sequence_);
  auto pam_sites = pam::search_pam_sites_fast(pam_seq_, genome_bits, params);

  // No PAM occurrences at all: skip tree construction silently.
  // (Some chromosomes — e.g. decoy sequences — may genuinely contain no
  // sites; treat this as a no-op rather than an error so batch indexing
  // does not abort on a single empty contig.)
  if (pam_sites.empty())
    return;

  // --- 2. Sequence extraction ------------------------------------------
  extract_sequences(pam_sites);

  if (leaves_.empty())
    throw std::runtime_error(
        "All PAM sites were discarded (contained N) for chromosome '" +
        chr_name_ + "'");

  // --- 3. Lexicographic sort -------------------------------------------
  std::sort(leaves_.begin(), leaves_.end(),
            [](const TSTLeaf &a, const TSTLeaf &b) {
              return a.guide_seq < b.guide_seq;
            });

  // --- 4 & 5. Partition → insert → serialize ---------------------------
  save();
}

// =========================================================================
// Sequence extraction
// =========================================================================

std::vector<uint8_t>
TernarySearchTree::encode_pam_bytes(std::string_view pam_str) {
  const int n = static_cast<int>(pam_str.size());
  std::vector<uint8_t> out((n + 1) / 2, 0);

  for (int i = 0; i < n; ++i) {
    uint8_t enc = iupac::encode_genome(pam_str[i]);
    if (i % 2 == 0)
      out[i / 2] = static_cast<uint8_t>(enc << 4); // high nibble
    else
      out[i / 2] |= enc; // low  nibble
  }
  return out;
}

void TernarySearchTree::extract_forward(int pos,
                                        std::vector<TSTLeaf> &dest) const {
  const int window = pam_length_ + max_bulges_;

  if (pos < 0 || pos + window > static_cast<int>(sequence_.size()))
    return;

  std::string_view window_view(sequence_.data() + pos, window);

  if (window_view.find('N') != std::string_view::npos)
    return;

  TSTLeaf leaf;
  leaf.guide_index = pos; // positive → forward strand

  if (!pam_at_start_) {
    // PAM at 3' end: window = [guide][pam+bulge_extra]
    // The guide is reversed before insertion (TST search order is 3'→5').
    std::string guide_raw(window_view.substr(0, guide_length_ + max_bulges_));
    std::reverse(guide_raw.begin(), guide_raw.end());
    leaf.guide_seq = std::move(guide_raw);

    // PAM: rightmost pam_limit_ chars, reversed for storage.
    std::string pam_raw(window_view.substr(guide_length_, pam_limit_));
    std::reverse(pam_raw.begin(), pam_raw.end());
    leaf.pam_seq_enc = encode_pam_bytes(pam_raw);
  } else {
    // PAM at 5' start: window = [pam][guide+bulge_extra]
    std::string pam_raw(window_view.substr(0, pam_limit_));
    std::reverse(pam_raw.begin(), pam_raw.end());
    leaf.pam_seq_enc = encode_pam_bytes(pam_raw);

    leaf.guide_seq = std::string(
        window_view.substr(pam_limit_, guide_length_ + max_bulges_));
  }

  dest.push_back(std::move(leaf));
}

void TernarySearchTree::extract_reverse(int pos,
                                        std::vector<TSTLeaf> &dest) const {
  const int window = pam_length_ + max_bulges_;

  if (pos < 0 || pos + window > static_cast<int>(sequence_.size()))
    return;

  std::string_view window_view(sequence_.data() + pos, window);

  if (window_view.find('N') != std::string_view::npos)
    return;

  std::string rc = reverse_complement(window_view);

  TSTLeaf leaf;
  leaf.guide_index = -pos; // negative → reverse strand

  if (!pam_at_start_) {
    // After RC, layout is the same as a forward hit (guide then pam).
    // For reverse strand + PAM-at-end the legacy code does NOT re-reverse
    // the guide; the RC already inverted orientation.
    std::string guide_raw(rc.substr(0, guide_length_ + max_bulges_));
    leaf.guide_seq = std::move(guide_raw);

    std::string pam_raw(rc.substr(guide_length_, pam_limit_));
    leaf.pam_seq_enc = encode_pam_bytes(pam_raw);
  } else {
    // PAM-at-start, reverse strand.
    std::string pam_raw(rc.substr(0, pam_limit_));
    std::reverse(pam_raw.begin(), pam_raw.end());
    leaf.pam_seq_enc = encode_pam_bytes(pam_raw);

    leaf.guide_seq =
        std::string(rc.substr(pam_limit_, guide_length_ + max_bulges_));
  }

  dest.push_back(std::move(leaf));
}

void TernarySearchTree::extract_sequences(const std::vector<int> &pam_sites) {
  leaves_.reserve(pam_sites.size());

  for (int site : pam_sites) {
    if (pam_at_start_) {
      // PAM-at-start convention (e.g. Cas12a):
      //   negative site → forward strand
      //   positive site → reverse strand
      if (site < 0)
        extract_forward(-site, leaves_);
      else
        extract_reverse(site, leaves_);
    } else {
      // PAM-at-end convention (e.g. SpCas9):
      //   positive site → forward strand
      //   negative site → reverse strand
      if (site > 0)
        extract_forward(site, leaves_);
      else
        extract_reverse(-site, leaves_);
    }
  }

  leaves_.shrink_to_fit();
}

// =========================================================================
// TST insertion
// =========================================================================

int TernarySearchTree::alloc_node() {
  int idx = nodes_used_++;
  if (idx >= static_cast<int>(nodes_.size()))
    nodes_.emplace_back();
  return idx;
}

void TernarySearchTree::insert(std::string_view guide_str, int leaf_idx,
                               int chunk_offset) {
  assert(!guide_str.empty());
  assert(nodes_used_ > 0 && "root must be allocated before first insert");

  const char *s = guide_str.data();
  const int encoded_leaf = -((leaf_idx - chunk_offset) + 1);

  int cur = 0; // start at root

  while (nodes_used_ > 0) {
    TSTNode &node = nodes_[cur];
    int d = static_cast<int>(static_cast<unsigned char>(*s)) -
            static_cast<int>(static_cast<unsigned char>(node.splitchar));

    if (d == 0) {
      ++s;
      if (*s == '\0') {
        leaves_[leaf_idx].next = node.eqkid;
        node.eqkid = encoded_leaf;
        return;
      }
      if (node.eqkid == 0) {
        node.eqkid = alloc_node();
        break;
      }
      cur = node.eqkid;
    } else if (d < 0) {
      if (node.lokid == 0) {
        node.lokid = alloc_node();
        break;
      }
      cur = node.lokid;
    } else {
      if (node.hikid == 0) {
        node.hikid = alloc_node();
        break;
      }
      cur = node.hikid;
    }
  }

  // Append new nodes for the remaining characters.
  while (true) {
    TSTNode &node = nodes_[cur];
    node.splitchar = *s;
    node.splitchar_enc = iupac::encode_genome(*s);

    ++s;
    if (*s == '\0') {
      node.eqkid = encoded_leaf;
      return;
    }
    node.eqkid = alloc_node();
    cur = node.eqkid;
  }
}

void TernarySearchTree::insert_balanced(int lo, int hi, int chunk_offset) {
  if (hi < lo)
    return;

  int mid = lo + (hi - lo) / 2;
  insert(leaves_[mid].guide_seq, mid, chunk_offset);
  insert_balanced(lo, mid - 1, chunk_offset);
  insert_balanced(mid + 1, hi, chunk_offset);
}

// ===========================================================================
// Serialization
// ===========================================================================

static uint8_t char_to_node_nibble(char c) {
  if (c == '0')
    return NULL_CHILD_NIBBLE;
  if (c == '_')
    return SENTINEL_NIBBLE;
  return iupac::encode_genome(c);
}

void TernarySearchTree::flush_pair(const char buf[2], int &buf_pos,
                                   std::ofstream &out) {
  uint8_t high = char_to_node_nibble(buf[0]);
  uint8_t low = char_to_node_nibble(buf[1]);

  uint8_t byte;
  if (buf[0] == '_')
    byte = static_cast<uint8_t>((SENTINEL_NIBBLE << 4) | low_nibble(low));
  else
    byte = pack_nibbles(high, low);

  out.put(static_cast<char>(byte));
  buf_pos = 0;
}

void TernarySearchTree::buffer_char(char c, char buf[2], int &buf_pos,
                                    std::ofstream &out) {
  buf[buf_pos++] = c;
  if (buf_pos == 2)
    flush_pair(buf, buf_pos, out);
}

void TernarySearchTree::serialize_node(int node_idx, std::ofstream &out,
                                       char buf[2], int &buf_pos) const {
  const TSTNode &node = nodes_[node_idx];

  buffer_char(node.splitchar, buf, buf_pos, out);

  if (node.lokid > 0)
    serialize_node(node.lokid, out, buf, buf_pos);
  else
    buffer_char('0', buf, buf_pos, out);

  if (node.hikid > 0)
    serialize_node(node.hikid, out, buf, buf_pos);
  else
    buffer_char('0', buf, buf_pos, out);

  if (node.eqkid > 0) {
    serialize_node(node.eqkid, out, buf, buf_pos);
  } else {
    buffer_char('_', buf, buf_pos, out);
    // Force flush so the 4-byte leaf pointer starts on a clean byte.
    if (buf_pos == 1) {
      buf[1] = '0';
      buf_pos = 2;
      flush_pair(buf, buf_pos, out);
    }
    out.write(reinterpret_cast<const char *>(&node.eqkid), sizeof(int));
  }
}

void TernarySearchTree::write_partition(int part, int chunk_start,
                                        int chunk_end) const {
  const int chunk_size = chunk_end - chunk_start;

  // Filename: <outdir_>/<pam_seq_>_<chr_name_>_<part>.bin
  const std::string filename = outdir_ + "/" + pam_seq_ + "_" + chr_name_ +
                               "_" + std::to_string(part) + ".bin";

  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out.is_open())
    throw std::runtime_error("Cannot open output file: " + filename);

  // ---- header ----
  out.write(reinterpret_cast<const char *>(&chunk_size), sizeof(int));
  out.write(reinterpret_cast<const char *>(&guide_length_), sizeof(int));

  // ---- leaf array ----
  for (int i = chunk_start; i < chunk_end; ++i) {
    const TSTLeaf &leaf = leaves_[i];

    out.write(reinterpret_cast<const char *>(&leaf.guide_index), sizeof(int));
    out.write(reinterpret_cast<const char *>(leaf.pam_seq_enc.data()),
              static_cast<std::streamsize>(leaf.pam_seq_enc.size()));

    if (leaf.next == 0) {
      out.put('0');
    } else {
      out.put('_');
      out.write(reinterpret_cast<const char *>(&leaf.next), sizeof(int));
    }
  }

  // ---- node count ----
  out.write(reinterpret_cast<const char *>(&nodes_used_), sizeof(int));

  // ---- serialized TST ----
  char buf[2] = {'\0', '\0'};
  int buf_pos = 0;
  serialize_node(0, out, buf, buf_pos);

  if (buf_pos == 1) {
    buf[1] = '0';
    buf_pos = 2;
    flush_pair(buf, buf_pos, out);
  }

  out.close();
  std::cout << "Written: " << filename << " (" << chunk_size << " leaves, "
            << nodes_used_ << " nodes)\n";
}

void TernarySearchTree::save() const {
  const int total = static_cast<int>(leaves_.size());
  const int groups = static_cast<int>(
      std::ceil(static_cast<double>(total) / LEAVES_PER_GROUP));

  for (int g = 0; g < groups; ++g) {
    const int chunk_start = g * LEAVES_PER_GROUP;
    const int chunk_end = std::min((g + 1) * LEAVES_PER_GROUP, total);

    // Reset node pool for this partition (safe const_cast: only mutable
    // build-state is touched, not the observable leaf data).
    auto *mut = const_cast<TernarySearchTree *>(this);
    mut->nodes_.clear();
    mut->nodes_.resize(static_cast<size_t>(chunk_end - chunk_start) *
                       pam_length_);
    mut->nodes_used_ = 0;

    mut->alloc_node(); // allocate root at index 0
    mut->nodes_[0].splitchar = '\0';
    mut->nodes_[0].splitchar_enc = 0;

    mut->insert_balanced(chunk_start, chunk_end - 1, chunk_start);

    write_partition(g + 1, chunk_start, chunk_end);
  }
}

// ===========================================================================
// Free function — pybind11 entry point
// ===========================================================================

void build_tree(const std::string &sequence, const std::string &chr_name,
                const std::string &pam_seq, int pam_length, int pam_limit,
                bool pam_at_start, const std::string &outdir, int max_bulges,
                int num_threads) {
  TernarySearchTree tst(sequence, chr_name, pam_seq, pam_length, pam_limit,
                        pam_at_start, outdir, max_bulges, num_threads);
  tst.build();
}

} // namespace crispritz
