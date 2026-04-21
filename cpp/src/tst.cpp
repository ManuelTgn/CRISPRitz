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

namespace crispritz
{

    // ===========================================================================
    // Constructor
    // ===========================================================================

    TernarySearchTree::TernarySearchTree(std::string_view sequence, std::string_view chr_name,
                                         std::string_view pam_seq, int pam_length, int pam_limit,
                                         bool pam_at_start, int max_bulges, int num_threads)
        : sequence_(sequence), chr_name_(chr_name), pam_seq_(pam_seq), pam_length_(pam_length),
          pam_limit_(pam_limit), guide_length_(pam_length - pam_limit), pam_at_start_(pam_at_start),
          max_bulges_(max_bulges), num_threads_(num_threads > 0 ? num_threads : 1)
    {
        if (guide_length_ <= 0)
            throw std::runtime_error("guide_length must be positive");
        if (pam_limit_ <= 0)
            throw std::runtime_error("pam_limit must be positive");
    }

    // ===========================================================================
    // Public entry point
    // ===========================================================================

    void TernarySearchTree::build()
    {
        // search PAM occurrences in input sequence
        pam::SearchParams params(pam_length_, pam_limit_, pam_at_start_, num_threads_);
        pam::CompactGenome genome_bits(sequence_);
        // vector<int> (PAM positions +/- strands)
        auto pam_sites = pam::search_pam_sites_fast(pam_seq_, genome_bits, params);

        // unlikely to fallback here -> no PAM occurrence found in input sequence?
        if (pam_sites.empty())
            return // if condition hit, skip tree construction

                // extract sequences from genome starting at PAM occurrence positions
                extract_sequences(pam_sites);

        // unlikely to fallback here -> all extracted sequences contain 'N'?
        // TODO: check if use cases triggering this behavior exist
        if (leaves_.empty())
            throw std::runtime_error("All PAM sites were discarded (contained N) for chromosome '" +
                                     chr_name_ + "'");

        // sort guide sequences lexicographically
        std::sort(leaves_.begin(), leaves_.end(),
                  [](const TSTLeaf& a, const TSTLeaf& b) { return a.guide_seq < b.guide_seq; });

        // Partition -> insert -> serialize
        // The legacy code splits large chromosomes into chunks of LEAVES_PER_GROUP
        // so that individual .bin files stay manageable.  Each chunk gets its own
        // TST built from scratch; the node pool is reset between chunks.
        save();
    }

    // ===========================================================================
    // Sequence extraction
    // ===========================================================================

    /**
     * Encode a PAM character substring into bit-packed bytes.
     *
     * Two IUPAC 4-bit codes per byte, high nibble first.  If the PAM length is
     * odd, the last character occupies only the high nibble of the last byte and
     * the low nibble is zero.
     *
     * Legacy reference: the inner do-while loop in saveTST() in mainParallel.cpp
     * that writes targetOnDNA[i].pamDNA.
     */
    std::vector<uint8_t> TernarySearchTree::encode_pam_bytes(std::string_view pam_str)
    {
        const int n = static_cast<int>(pam_str.size());
        // ceil(n / 2) bytes: each byte holds two nibbles
        std::vector<uint8_t> out((n + 1) / 2, 0);

        for (int i = 0; i < n; ++i)
        {
            uint8_t enc = iupac::encode_genome(pam_str[i]);
            if (i % 2 == 0)
                out[i / 2] = static_cast<uint8_t>(enc << 4); // high nibble
            else
                out[i / 2] |= enc; // low nibble
        }
        return out;
    }

    void TernarySearchTree::extract_forward(int pos, std::vector<TSTLeaf>& dest) const
    {
        // how many characters to read: guide + PAM + extra bases for bulges
        const int window = pam_length_ + max_bulges_;

        // bounds check: site must leave enough room for the full window
        if (pos < 0 || pos + window > static_cast<int>(sequence_.size()))
            return;

        std::string_view window_view(sequence_.data() + pos, window);

        // discard any site that contains 'N'
        if (window_view.find('N') != std::string_view::npos)
            return;

        TSTLeaf leaf;
        leaf.guide_index = pos; // positive = forward strand

        if (!pam_at_start_)
        {
            // PAM is at the 3' end (right side of the window).
            // Legacy layout: window = [guide(0..guide_length_-1)][pam(guide_length_..)]
            //
            // The guide substring is REVERSED before insertion so that the TST
            // nodes are ordered 3'->5' (matching the search direction in
            // nearsearch()).
            //
            // The PAM substring is also reversed, then stored; the saved order
            // will be reversed again when building the PAM string during search.
            std::string guide_raw(window_view.substr(0, guide_length_ + max_bulges_));
            std::reverse(guide_raw.begin(), guide_raw.end());
            leaf.guide_seq = std::move(guide_raw);

            // PAM: the rightmost pam_limit_ chars of the original window,
            // stored as-is (the reversal happened to the guide only).
            std::string pam_raw(window_view.substr(guide_length_, pam_limit_));
            std::reverse(pam_raw.begin(), pam_raw.end());
            leaf.pam_seq_enc = encode_pam_bytes(pam_raw);
        }
        else
        {
            // PAM is at the 5' start (left side of the window).
            // Legacy layout: window = [pam(0..pam_limit_-1)][guide(pam_limit_..)]
            //
            // The guide substring is taken as-is (already 5'->3').
            // The PAM is reversed for storage (legacy reversal in mainParallel.cpp
            // lines: tmp_pam_str = target.substr(0, pamRNA.length()); reverse(...))
            std::string pam_raw(window_view.substr(0, pam_limit_));
            std::reverse(pam_raw.begin(), pam_raw.end());
            leaf.pam_seq_enc = encode_pam_bytes(pam_raw);

            leaf.guide_seq =
                std::string(window_view.substr(pam_limit_, guide_length_ + max_bulges_));
        }

        dest.push_back(std::move(leaf));
    }

    void TernarySearchTree::extract_reverse(int pos, std::vector<TSTLeaf>& dest) const
    {
        // pos is the raw absolute position (the PAM search returns it negated for
        // reverse-strand hits; the caller passes the absolute value here).
        const int window = pam_length_ + max_bulges_;

        if (pos < 0 || pos + window > static_cast<int>(sequence_.size()))
            return;

        std::string_view window_view(sequence_.data() + pos, window);

        if (window_view.find('N') != std::string_view::npos)
            return;

        // reverse complement the entire window first
        std::string rc = reverse_complement(window_view);

        TSTLeaf leaf;
        leaf.guide_index = -pos; // negative = reverse strand (legacy convention)

        if (!pam_at_start_)
        {
            // after RC, the window is now [guide][pam] as seen on the opposite
            // strand.  Apply the same reversal of the guide as for forward
            std::string guide_raw(rc.substr(0, guide_length_ + max_bulges_));
            // for the reverse strand with PAM-at-end, the legacy code does NOT
            // reverse again (the RC already inverted orientation).
            leaf.guide_seq = std::move(guide_raw);

            std::string pam_raw(rc.substr(guide_length_, pam_limit_));
            leaf.pam_seq_enc = encode_pam_bytes(pam_raw);
        }
        else
        {
            // PAM-at-start, reverse strand:
            // after RC, pam occupies the first pam_limit_ characters.
            std::string pam_raw(rc.substr(0, pam_limit_));
            std::reverse(pam_raw.begin(), pam_raw.end());
            leaf.pam_seq_enc = encode_pam_bytes(pam_raw);

            leaf.guide_seq = std::string(rc.substr(pam_limit_, guide_length_ + max_bulges_));
        }

        dest.push_back(std::move(leaf));
    }

    void TernarySearchTree::extract_sequences(const std::vector<int>& pam_sites)
    {
        // pre-allocate generously; many sites will be discarded due to Ns
        leaves_.reserve(pam_sites.size());

        for (int site : pam_sites)
        {
            if (pam_at_start_)
            {
                // positive index -> PAM-at-start convention: negative site =
                // forward strand, positive site = reverse strand
                // (Legacy mainParallel.cpp: if (pamIndices[i] < 0) -> positive strand)
                if (site < 0)
                    extract_forward(-site, leaves_);
                else
                    extract_reverse(site, leaves_);
            }
            else
            {
                // PAM-at-end convention: positive site = forward strand,
                // negative site = reverse strand
                if (site > 0)
                    extract_forward(site, leaves_);
                else
                    extract_reverse(-site, leaves_);
            }
        }

        leaves_.shrink_to_fit();
    }

    // ===========================================================================
    // TST insertion
    // ===========================================================================

    int TernarySearchTree::alloc_node()
    {
        int idx = nodes_used_++;
        if (idx >= static_cast<int>(nodes_.size()))
            nodes_.emplace_back(); // zero-initialised by TSTNode default ctor
        return idx;
    }

    void TernarySearchTree::insert(std::string_view guide_str, int leaf_idx, int chunk_offset)
    {
        assert(!guide_str.empty());
        assert(nodes_used_ > 0 && "root must be allocated before first insert");

        const char* s = guide_str.data();

        // within-chunk leaf index encoded in the eqkid field (legacy: (i2+1)*-1)
        const int encoded_leaf = -((leaf_idx - chunk_offset) + 1);

        int cur = 0; // start from root

        while (nodes_used_ > 1) // at least root exists
        {
            TSTNode& node = nodes_[cur];
            int d = static_cast<int>(static_cast<unsigned char>(*s)) -
                    static_cast<int>(static_cast<unsigned char>(node.splitchar));

            if (d == 0)
            {
                ++s;
                if (*s == '\0')
                {
                    // end of string: link collision chain
                    leaves_[leaf_idx].next = node.eqkid;
                    node.eqkid = encoded_leaf;
                    return;
                }
                if (node.eqkid == 0)
                {
                    node.eqkid = alloc_node();
                    break;
                }
                cur = node.eqkid;
            }
            else if (d < 0)
            {
                if (node.lokid == 0)
                {
                    node.lokid = alloc_node();
                    break;
                }
                cur = node.lokid;
            }
            else
            {
                if (node.hikid == 0)
                {
                    node.hikid = alloc_node();
                    break;
                }
                cur = node.hikid;
            }
        }

        // append new node(s) for the remaining characters of the string
        while (true)
        {
            TSTNode& node = nodes_[cur];
            node.splitchar = *s;
            node.splitchar_enc = iupac::encode_genome(*s);

            ++s;
            if (*s == '\0')
            {
                node.eqkid = encoded_leaf;
                return;
            }
            node.eqkid = alloc_node();
            cur = node.eqkid;
        }
    }

    void TernarySearchTree::insert_balanced(int lo, int hi, int chunk_offset)
    {
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

    /**
     * Encode an IUPAC character to its 4-bit code for the node stream.
     *
     * '0' (null child) maps to 0x00.
     * '_' (end-of-node sentinel) maps to 0b1111 (the high nibble 0b1111 triggers
     * the sentinel detection in the legacy readPair / deSerialize code).
     *
     * All other characters use iupac::encode_genome.
     */
    static uint8_t char_to_node_nibble(char c)
    {
        if (c == '0')
            return NULL_CHILD_NIBBLE;
        if (c == '_')
            return SENTINEL_NIBBLE;
        return iupac::encode_genome(c);
    }

    void TernarySearchTree::flush_pair(const char buf[2], int& buf_pos, std::ofstream& out)
    {
        uint8_t high = char_to_node_nibble(buf[0]);
        uint8_t low = char_to_node_nibble(buf[1]);

        // special case: the sentinel '_' in the high nibble is written as 0b1111
        // (legacy writePair: if (pairNuc[0] == '_') bitNuc = 0b1111).
        uint8_t byte;
        if (buf[0] == '_')
            byte = static_cast<uint8_t>((SENTINEL_NIBBLE << 4) | low_nibble(low));
        else
            byte = pack_nibbles(high, low);

        out.put(static_cast<char>(byte));
        buf_pos = 0;
    }

    void TernarySearchTree::buffer_char(char c, char buf[2], int& buf_pos, std::ofstream& out)
    {
        buf[buf_pos++] = c;
        if (buf_pos == 2)
            flush_pair(buf, buf_pos, out);
    }

    void TernarySearchTree::serialize_node(int node_idx, std::ofstream& out, char buf[2],
                                           int& buf_pos) const
    {
        const TSTNode& node = nodes_[node_idx];

        // write this node's character into the buffer
        buffer_char(node.splitchar, buf, buf_pos, out);

        // lokid subtree (or null sentinel)
        if (node.lokid > 0)
            serialize_node(node.lokid, out, buf, buf_pos);
        else
            buffer_char('0', buf, buf_pos, out);

        // hikid subtree (or null sentinel)
        if (node.hikid > 0)
            serialize_node(node.hikid, out, buf, buf_pos);
        else
            buffer_char('0', buf, buf_pos, out);

        // eqkid: either a subtree or a leaf pointer
        if (node.eqkid > 0)
        {
            serialize_node(node.eqkid, out, buf, buf_pos);
        }
        else
        {
            // leaf pointer: write sentinel '_' then flush, then the raw 4-byte int
            buffer_char('_', buf, buf_pos, out);
            // if '_' was the first nibble, the sentinel byte may not have been
            // flushed yet.  Force a flush with a dummy second nibble so the int
            // that follows starts on a clean byte boundary.
            if (buf_pos == 1)
            {
                buf[1] = '0'; // dummy low nibble
                buf_pos = 2;
                flush_pair(buf, buf_pos, out);
            }
            // write the signed leaf-pointer int (4 bytes, platform endianness -
            // matching legacy fileTree.write((char*)&p->eqkid, sizeof(int)))
            out.write(reinterpret_cast<const char*>(&node.eqkid), sizeof(int));
        }
    }

    void TernarySearchTree::write_partition(int part, int chunk_start, int chunk_end) const
    {
        const int chunk_size = chunk_end - chunk_start;

        // filename: <pam_seq>_<chr_name>_<part>.bin  (legacy naming convention)
        std::string filename = pam_seq_ + "_" + chr_name_ + "_" + std::to_string(part) + ".bin";
        std::ofstream out(filename, std::ios::out | std::ios::binary);
        if (!out.is_open())
            throw std::runtime_error("Cannot open output file: " + filename);

        // ---- header ----
        out.write(reinterpret_cast<const char*>(&chunk_size), sizeof(int));
        out.write(reinterpret_cast<const char*>(&guide_length_), sizeof(int));

        // ---- leaf array ----
        for (int i = chunk_start; i < chunk_end; ++i)
        {
            const TSTLeaf& leaf = leaves_[i];

            // guide genomic index (signed 4-byte int)
            out.write(reinterpret_cast<const char*>(&leaf.guide_index), sizeof(int));

            // bit-packed PAM bytes
            out.write(reinterpret_cast<const char*>(leaf.pam_seq_enc.data()),
                      static_cast<std::streamsize>(leaf.pam_seq_enc.size()));

            // next-leaf link
            if (leaf.next == 0)
            {
                out.put('0');
            }
            else
            {
                out.put('_');
                out.write(reinterpret_cast<const char*>(&leaf.next), sizeof(int));
            }
        }

        // ---- node count ----
        out.write(reinterpret_cast<const char*>(&nodes_used_), sizeof(int));

        // ---- serialized TST ----
        char buf[2] = {'\0', '\0'};
        int buf_pos = 0;
        serialize_node(0, out, buf, buf_pos);

        // if an odd number of characters was written, the last nibble is still
        // buffered.  Flush it with a null second nibble
        if (buf_pos == 1)
        {
            buf[1] = '0';
            buf_pos = 2;
            flush_pair(buf, buf_pos, out);
        }

        out.close();
        std::cout << "Written: " << filename << " (" << chunk_size << " leaves, " << nodes_used_
                  << " nodes)\n";
    }

    void TernarySearchTree::save() const
    {
        const int total = static_cast<int>(leaves_.size());
        const int groups =
            static_cast<int>(std::ceil(static_cast<double>(total) / LEAVES_PER_GROUP));

        for (int g = 0; g < groups; ++g)
        {
            const int chunk_start = g * LEAVES_PER_GROUP;
            const int chunk_end = std::min((g + 1) * LEAVES_PER_GROUP, total);

            // ---- Reset node pool for this partition ----
            // We cast away const here only on the mutable node state; leaves_ is
            // not modified.  A cleaner design would separate build-state from the
            // const-observable interface, but casting is safe here because save()
            // is always called from build() on a fully constructed object.
            auto* mutable_this = const_cast<TernarySearchTree*>(this);
            mutable_this->nodes_.clear();
            mutable_this->nodes_.resize(static_cast<size_t>(chunk_end - chunk_start) * pam_length_);
            mutable_this->nodes_used_ = 0;

            // allocate root
            mutable_this->alloc_node();
            nodes_[0].splitchar = '\0';
            nodes_[0].splitchar_enc = 0;

            mutable_this->insert_balanced(chunk_start, chunk_end - 1, chunk_start);

            write_partition(g + 1, chunk_start, chunk_end);
        }
    }

    // ===========================================================================
    // Free function (pybind11 entry point)
    // ===========================================================================

    void build_tree(const std::string& sequence, const std::string& chr_name,
                    const std::string& pam_seq, int pam_length, int pam_limit, bool pam_at_start,
                    int max_bulges, int num_threads)
    {
        TernarySearchTree tst(sequence, chr_name, pam_seq, pam_length, pam_limit, pam_at_start,
                              max_bulges, num_threads);
        tst.build();
    }

} // namespace crispritz
