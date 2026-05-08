#include "search.hpp"

#include "offtarget.hpp"
#include "tst_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace crispritz
{

    // =========================================================================
    // Internal helpers
    // =========================================================================

    namespace
    {

        // ---------------------------------------------------------------------
        // Decode a packed-nibble byte back to a printable IUPAC character.
        // The lookup table is the same as convertBitsetToChar() in the legacy 
        // code.
        // ---------------------------------------------------------------------
        static const char NUC_DECODE[16] = {
            '-', // 0000  gap / bulge placeholder
            'A', // 0001
            'C', // 0010
            '?', // 0011  M – should not appear in mismatch context
            'G', // 0100
            'R', // 0101
            'S', // 0110
            'V', // 0111
            'T', // 1000
            'W', // 1001
            'Y', // 1010
            'H', // 1011
            'K', // 1100
            'D', // 1101
            'B', // 1110
            'N', // 1111
        };

        /** @brief True iff a & b share at least one set bit (IUPAC match). */
        inline bool nuc_match(uint8_t a, uint8_t b) noexcept
        {
            return (a & b) != 0;
        }

        /** @brief True iff both nibbles are non-zero and share NO set bit (mismatch). */
        inline bool nuc_mismatch(uint8_t a, uint8_t b) noexcept
        {
            return (a != 0) && (b != 0) && ((a & b) == 0);
        }

        // ---------------------------------------------------------------------------
        // Build a printable PAM string from a leaf's packed PAM bytes.
        // The pam_limit first nibbles (high-nibble first within each byte) are decoded.
        // ---------------------------------------------------------------------------
        static std::string decode_pam(const std::vector<uint8_t>& pam_enc, int pam_limit)
        {
            std::string result(pam_limit, 'N');
            for (int i = 0; i < pam_limit; ++i)
            {
                uint8_t nibble =
                    (i % 2 == 0) ? high_nibble(pam_enc[i / 2]) : low_nibble(pam_enc[i / 2]);
                result[i] = NUC_DECODE[nibble & 0x0F];
            }
            return result;
        }

        // ---------------------------------------------------------------------------
        // Per-thread working state for the recursive TST search.
        // Each thread maintains its own instance to avoid false sharing.
        // ---------------------------------------------------------------------------
        struct ThreadState
        {
            // Current guide character index (position in the guide we are matching)
            int gi = 0;
            // Current target character index (position we have written into target_buf)
            int ti = 0;

            // Buffers accumulating the guide / target strings as we descend the TST.
            // Size = guide_len + max_dna_bulges + 1 (extra +1 for null sentinel).
            std::vector<char> guide_buf;
            std::vector<uint8_t> guide_enc; // 4-bit encoding of each guide_buf char
            std::vector<char> target_buf;
            std::vector<uint8_t> target_enc;

            // Collected hits for this thread (flushed into the shared result).
            std::vector<OffTarget> hits;

            void init(int guide_len, int max_dna)
            {
                int cap = guide_len + max_dna + 2;
                guide_buf.assign(cap, '\0');
                guide_enc.assign(cap, 0);
                target_buf.assign(cap, '\0');
                target_enc.assign(cap, 0);
                gi = ti = 0;
                hits.clear();
            }
        };

        // ---------------------------------------------------------------------------
        // Forward-declare the recursive near-search so the compiler sees it.
        // ---------------------------------------------------------------------------
        static void nearsearch(const TSTIndex& idx, int node_idx,
                               const uint8_t* guide_enc_full, // full encoded guide
                               int guide_pos,                 // position in guide
                               int mm_left, int dna_left, int rna_left,
                               int bulge_type_rna, // >0 if RNA bulge occurred
                               int bulge_type_dna, // >0 if DNA bulge occurred
                               bool go_lo_hi, const SearchParams& params,
                               const std::string& chr_name, int guide_idx,
                               const std::string& guide_str_original, int offset_guide_len,
                               ThreadState& ts);

        // ---------------------------------------------------------------------------
        // Called when the search has consumed the full guide length (ti == guide_end).
        // Visits every leaf reachable from node_idx (traversing lo/eq/hi children)
        // and records an OffTarget for each PAM site found.
        // ---------------------------------------------------------------------------
        static void collect_leaves(const TSTIndex& idx, int node_idx, int mm_used, int dna_used,
                                   int rna_used, int bulge_type_rna, int bulge_type_dna,
                                   const SearchParams& params, const std::string& chr_name,
                                   int guide_idx, const std::string& guide_str_original,
                                   int offset_guide_len, ThreadState& ts)
        {
            if (node_idx <= 0 || node_idx >= static_cast<int>(idx.nodes.size()))
                return;

            const TSTIndexNode& node = idx.nodes[node_idx];

            // Recurse into lo/hi children first (they may also lead to leaves).
            if (node.lokid > 0)
                collect_leaves(idx, node.lokid, mm_used, dna_used, rna_used, bulge_type_rna,
                               bulge_type_dna, params, chr_name, guide_idx, guide_str_original,
                               offset_guide_len, ts);
            if (node.hikid > 0)
                collect_leaves(idx, node.hikid, mm_used, dna_used, rna_used, bulge_type_rna,
                               bulge_type_dna, params, chr_name, guide_idx, guide_str_original,
                               offset_guide_len, ts);

            if (node.eqkid >= 0)
            {
                // Internal node – follow equal child.
                if (node.eqkid > 0)
                    collect_leaves(idx, node.eqkid, mm_used, dna_used, rna_used, bulge_type_rna,
                                   bulge_type_dna, params, chr_name, guide_idx, guide_str_original,
                                   offset_guide_len, ts);
                return;
            }

            // ---- This node IS a leaf ----
            int leaf_chain = (-node.eqkid - 1);

            // Build the guide string (with PAM Ns) in canonical display orientation.
            const int guide_len = params.guide_length();
            const int dna_total = params.max_dna_bulges;

            // The guide_buf holds the guide characters in TST traversal order.
            // For PAM-at-end (pam_at_start==false) the guide was reversed on insertion,
            // so we reverse it back for display.
            std::string g_display(ts.guide_buf.data(),
                                  guide_len + dna_total - (params.max_dna_bulges - dna_used));
            std::string t_display(ts.target_buf.data(),
                                  guide_len + dna_total - (params.max_dna_bulges - dna_used));

            if (!params.pam_at_start)
            {
                std::reverse(g_display.begin(), g_display.end());
                std::reverse(t_display.begin(), t_display.end());
            }

            // Walk the linked list of PAM sites for this leaf.
            while (leaf_chain >= 0 && leaf_chain < static_cast<int>(idx.leaves.size()))
            {
                const TSTIndexLeaf& leaf = idx.leaves[leaf_chain];

                // Attach PAM Ns to the guide display string.
                std::string g_full = g_display;
                std::string t_full = t_display;

                std::string pam_str = decode_pam(leaf.pam_enc, params.pam_limit);

                if (params.pam_at_start)
                {
                    g_full = std::string(params.pam_limit, 'N') + g_full;
                    t_full = pam_str + t_full;
                }
                else
                {
                    g_full += std::string(params.pam_limit, 'N');
                    t_full += pam_str;
                }

                // Determine bulge type label.
                std::string btype;
                if (bulge_type_rna == 0 && bulge_type_dna == 0)
                    btype = "X";
                else if (bulge_type_rna != 0 && bulge_type_dna != 0)
                    btype = "DNA,RNA";
                else if (bulge_type_dna != 0)
                    btype = "DNA";
                else
                    btype = "RNA";

                int bsize = dna_used + rna_used;

                // Resolve strand and genomic position.
                char strand;
                int gpos, cpos;

                if (leaf.guide_index < 0)
                {
                    // Negative index → stored as negative, strand depends on pam_at_start.
                    int abs_pos = -leaf.guide_index;
                    gpos = abs_pos;
                    cpos = abs_pos;
                    strand = params.pam_at_start ? '+' : '-';
                }
                else
                {
                    // Positive index → forward strand position in legacy convention.
                    int raw = leaf.guide_index + params.max_dna_bulges - dna_used + rna_used +
                              offset_guide_len;
                    gpos = raw;

                    if (bulge_type_rna == 0 && bulge_type_dna == 0)
                        cpos = raw;
                    else if (bulge_type_rna != 0 && bulge_type_dna == 0)
                        cpos = raw - rna_used;
                    else if (bulge_type_rna == 0 && bulge_type_dna != 0)
                        cpos = raw + dna_used;
                    else
                        cpos = raw + dna_used - rna_used;

                    strand = params.pam_at_start ? '-' : '+';
                }

                OffTarget ot;
                ot.bulge_type = btype;
                ot.guide_seq = g_full;
                ot.target_seq = t_full;
                ot.chromosome = chr_name;
                ot.genomic_pos = gpos;
                ot.cluster_pos = cpos;
                ot.strand = strand;
                ot.mismatches = mm_used;
                ot.bulge_size = bsize;
                ot.total_score = mm_used + bsize;

                ts.hits.push_back(std::move(ot));

                leaf_chain = leaf.next;
            }
        }

        // ---------------------------------------------------------------------------
        // Leaf-only variant of collect_leaves: called when the last character of
        // the guide has just been matched/mismatched and we are at a terminal node
        // (eqkid == 0, i.e. the node itself is the final node of a guide sequence).
        // ---------------------------------------------------------------------------
        static void collect_single_leaf(const TSTIndex& idx, int node_idx, int mm_used,
                                        int dna_used, int rna_used, int bulge_type_rna,
                                        int bulge_type_dna, const SearchParams& params,
                                        const std::string& chr_name, int guide_idx,
                                        const std::string& guide_str_original, int offset_guide_len,
                                        ThreadState& ts)
        {
            if (node_idx <= 0 || node_idx >= static_cast<int>(idx.nodes.size()))
                return;

            const TSTIndexNode& node = idx.nodes[node_idx];
            if (node.eqkid >= 0)
                return; // not a leaf

            int leaf_chain = (-node.eqkid - 1);

            // Same display logic as collect_leaves.
            const int guide_len = params.guide_length();
            const int dna_total = params.max_dna_bulges;
            int used_len = guide_len + dna_total - (params.max_dna_bulges - dna_used);

            std::string g_display(ts.guide_buf.data(), used_len);
            std::string t_display(ts.target_buf.data(), used_len);

            if (!params.pam_at_start)
            {
                std::reverse(g_display.begin(), g_display.end());
                std::reverse(t_display.begin(), t_display.end());
            }

            while (leaf_chain >= 0 && leaf_chain < static_cast<int>(idx.leaves.size()))
            {
                const TSTIndexLeaf& leaf = idx.leaves[leaf_chain];

                std::string g_full = g_display;
                std::string t_full = t_display;
                std::string pam_str = decode_pam(leaf.pam_enc, params.pam_limit);

                if (params.pam_at_start)
                {
                    g_full = std::string(params.pam_limit, 'N') + g_full;
                    t_full = pam_str + t_full;
                }
                else
                {
                    g_full += std::string(params.pam_limit, 'N');
                    t_full += pam_str;
                }

                std::string btype;
                if (bulge_type_rna == 0 && bulge_type_dna == 0)
                    btype = "X";
                else if (bulge_type_rna != 0 && bulge_type_dna != 0)
                    btype = "DNA,RNA";
                else if (bulge_type_dna != 0)
                    btype = "DNA";
                else
                    btype = "RNA";

                int bsize = dna_used + rna_used;

                char strand;
                int gpos, cpos;

                if (leaf.guide_index < 0)
                {
                    gpos = -leaf.guide_index;
                    cpos = gpos;
                    strand = params.pam_at_start ? '+' : '-';
                }
                else
                {
                    int raw = leaf.guide_index + params.max_dna_bulges - dna_used + rna_used +
                              offset_guide_len;
                    gpos = raw;

                    if (bulge_type_rna == 0 && bulge_type_dna == 0)
                        cpos = raw;
                    else if (bulge_type_rna != 0 && bulge_type_dna == 0)
                        cpos = raw - rna_used;
                    else if (bulge_type_rna == 0 && bulge_type_dna != 0)
                        cpos = raw + dna_used;
                    else
                        cpos = raw + dna_used - rna_used;

                    strand = params.pam_at_start ? '-' : '+';
                }

                OffTarget ot;
                ot.bulge_type = btype;
                ot.guide_seq = g_full;
                ot.target_seq = t_full;
                ot.chromosome = chr_name;
                ot.genomic_pos = gpos;
                ot.cluster_pos = cpos;
                ot.strand = strand;
                ot.mismatches = mm_used;
                ot.bulge_size = bsize;
                ot.total_score = mm_used + bsize;

                ts.hits.push_back(std::move(ot));

                leaf_chain = leaf.next;
            }
        }

        // ---------------------------------------------------------------------------
        // Core recursive near-search on a single TST partition.
        //
        // Parameters mirror nearsearch() in searchOnTST.cpp:
        //   node_idx    – current TST node index
        //   guide_pos   – position in the guide we are currently matching
        //   mm_left     – remaining mismatches allowed
        //   dna_left    – remaining DNA bulges allowed
        //   rna_left    – remaining RNA bulges allowed
        //   go_lo_hi    – whether to visit lo/hi children before eq
        // ---------------------------------------------------------------------------
        static void nearsearch(const TSTIndex& idx, int node_idx, const uint8_t* guide_enc_full,
                               int guide_pos, int mm_left, int dna_left, int rna_left,
                               int bulge_type_rna, int bulge_type_dna, bool go_lo_hi,
                               const SearchParams& params, const std::string& chr_name,
                               int guide_idx, const std::string& guide_str_original,
                               int offset_guide_len, ThreadState& ts)
        {
            if (node_idx <= 0 || node_idx >= static_cast<int>(idx.nodes.size()))
                return;

            const TSTIndexNode& node = idx.nodes[node_idx];

            // -- lo / hi traversal (same-depth siblings) --
            if (go_lo_hi)
            {
                if (node.lokid > 0)
                    nearsearch(idx, node.lokid, guide_enc_full, guide_pos, mm_left, dna_left,
                               rna_left, bulge_type_rna, bulge_type_dna, true, params, chr_name,
                               guide_idx, guide_str_original, offset_guide_len, ts);
                if (node.hikid > 0)
                    nearsearch(idx, node.hikid, guide_enc_full, guide_pos, mm_left, dna_left,
                               rna_left, bulge_type_rna, bulge_type_dna, true, params, chr_name,
                               guide_idx, guide_str_original, offset_guide_len, ts);
            }

            const int guide_len = params.guide_length();
            const int end_pos =
                guide_len + params.max_dna_bulges - (params.max_dna_bulges - dna_left);

            // -- Have we consumed the full guide length? → collect leaves --
            if (ts.ti == end_pos)
            {
                collect_leaves(idx, node_idx, params.max_mismatches - mm_left,
                               params.max_dna_bulges - dna_left, params.max_rna_bulges - rna_left,
                               bulge_type_rna, bulge_type_dna, params, chr_name, guide_idx,
                               guide_str_original, offset_guide_len, ts);
                return;
            }

            // -- Internal node: advance along equal child (consuming one guide char) --
            if (node.eqkid > 0)
            {
                uint8_t g_enc = guide_enc_full[guide_pos];
                char g_chr = static_cast<char>(
                    NUC_DECODE[g_enc & 0x0F] == '?' ? 'N' : NUC_DECODE[g_enc & 0x0F]);

                // Save current guide / target positions.
                ts.guide_buf[ts.gi] = g_chr;
                ts.guide_enc[ts.gi] = g_enc;
                ts.gi++;

                // ---- MATCH ----
                if (nuc_match(g_enc, node.splitchar_enc))
                {
                    ts.target_buf[ts.ti] = node.splitchar;
                    ts.target_enc[ts.ti] = node.splitchar_enc;
                    ts.ti++;
                    nearsearch(idx, node.eqkid, guide_enc_full, guide_pos + 1, mm_left, dna_left,
                               rna_left, bulge_type_rna, bulge_type_dna, true, params, chr_name,
                               guide_idx, guide_str_original, offset_guide_len, ts);
                    ts.ti--;
                }

                // ---- MISMATCH ----
                if (mm_left > 0)
                {
                    // Lower-case the target character to flag the mismatch.
                    char lower_t = static_cast<char>(node.splitchar | 0x20);
                    ts.target_buf[ts.ti] = lower_t;
                    ts.target_enc[ts.ti] = node.splitchar_enc;
                    ts.ti++;
                    nearsearch(idx, node.eqkid, guide_enc_full, guide_pos + 1, mm_left - 1,
                               dna_left, rna_left, bulge_type_rna, bulge_type_dna, true, params,
                               chr_name, guide_idx, guide_str_original, offset_guide_len, ts);
                    ts.ti--;
                }

                // ---- RNA BULGE (gap in target / extra base in guide) ----
                // The guide advances by 1 but the target gets a '-' gap character.
                if (rna_left > 0)
                {
                    ts.target_buf[ts.ti] = '-';
                    ts.target_enc[ts.ti] = 0x00; // gap encoding
                    ts.ti++;
                    nearsearch(idx, node_idx, guide_enc_full, guide_pos + 1, mm_left, dna_left,
                               rna_left - 1, bulge_type_rna + 1, bulge_type_dna, true, params,
                               chr_name, guide_idx, guide_str_original, offset_guide_len, ts);
                    ts.ti--;
                }

                // ---- DNA BULGE (gap in guide / extra base in target) ----
                // The target advances by 1 but the guide gets a '-' gap character.
                if (dna_left > 0)
                {
                    // Overwrite the guide character we just pushed with a gap.
                    ts.guide_buf[ts.gi - 1] = '-';
                    ts.guide_enc[ts.gi - 1] = 0x00;

                    ts.target_buf[ts.ti] = node.splitchar;
                    ts.target_enc[ts.ti] = node.splitchar_enc;
                    ts.ti++;
                    // guide_pos stays the same (we did NOT consume a guide character).
                    nearsearch(idx, node.eqkid, guide_enc_full, guide_pos, mm_left, dna_left - 1,
                               rna_left, bulge_type_rna, bulge_type_dna + 1, true, params, chr_name,
                               guide_idx, guide_str_original, offset_guide_len, ts);
                    ts.ti--;

                    // Restore guide character.
                    ts.guide_buf[ts.gi - 1] = g_chr;
                    ts.guide_enc[ts.gi - 1] = g_enc;
                }

                ts.gi--;
            }
            else if (node.eqkid < 0)
            {
                // ---- Terminal node: last character of the guide sequence ----
                // The node character must be matched (or mismatched within budget).
                if (guide_pos >= guide_len)
                    return; // guard: out-of-bounds guide_pos

                uint8_t g_enc = guide_enc_full[guide_pos];
                char g_chr = static_cast<char>(
                    NUC_DECODE[g_enc & 0x0F] == '?' ? 'N' : NUC_DECODE[g_enc & 0x0F]);

                ts.guide_buf[ts.gi] = g_chr;
                ts.guide_enc[ts.gi] = g_enc;
                ts.gi++;

                bool is_match = nuc_match(g_enc, node.splitchar_enc);
                char t_char = is_match ? node.splitchar : static_cast<char>(node.splitchar | 0x20);
                int mm_cost = is_match ? 0 : 1;

                if (mm_left - mm_cost >= 0)
                {
                    ts.target_buf[ts.ti] = t_char;
                    ts.target_enc[ts.ti] = node.splitchar_enc;
                    ts.ti++;

                    if (ts.ti == end_pos)
                    {
                        collect_single_leaf(
                            idx, node_idx, params.max_mismatches - (mm_left - mm_cost),
                            params.max_dna_bulges - dna_left, params.max_rna_bulges - rna_left,
                            bulge_type_rna, bulge_type_dna, params, chr_name, guide_idx,
                            guide_str_original, offset_guide_len, ts);
                    }

                    ts.ti--;
                }

                ts.gi--;
            }
        }

        // ---------------------------------------------------------------------------
        // Build the PositionalProfile for one guide from its collected OffTargets.
        // ---------------------------------------------------------------------------
        static PositionalProfile build_profile(const std::string& guide_str,
                                               const std::vector<OffTarget>& hits,
                                               const SearchParams& params)
        {
            PositionalProfile prof;
            prof.guide = guide_str + std::string(params.pam_limit, 'N');
            int guide_len = params.guide_length();
            int max_mm = params.max_mismatches;
            int max_dna = params.max_dna_bulges;
            int max_rna = params.max_rna_bulges;

            prof.mm_per_pos.assign(guide_len, 0);
            prof.dna_per_pos.assign(guide_len, 0);
            prof.rna_per_pos.assign(guide_len, 0);
            prof.mm_counts.assign(max_mm + 1, 0);

            prof.dna_counts.assign(max_mm + 1, std::vector<int>(max_dna > 0 ? max_dna : 1, 0));
            prof.rna_counts.assign(max_mm + 1, std::vector<int>(max_rna > 0 ? max_rna : 1, 0));

            prof.joint_counts.assign(
                max_mm + 1,
                std::vector<std::vector<int>>(max_dna + 1, std::vector<int>(max_rna + 1, 0)));

            prof.nuc_per_pos.assign(
                max_mm + 1, std::vector<std::vector<int>>(4, std::vector<int>(guide_len, 0)));

            for (const OffTarget& ot : hits)
            {
                int mm = ot.mismatches;
                int bd = 0, br = 0;

                if (ot.bulge_type == "DNA")
                    bd = ot.bulge_size;
                else if (ot.bulge_type == "RNA")
                    br = ot.bulge_size;
                else if (ot.bulge_type == "DNA,RNA")
                {
                    // Approximate split — we only have total in OffTarget;
                    // exact per-type counts require the detailed traversal arrays.
                    // Use the convention that bd + br == bulge_size.
                    bd = ot.bulge_size / 2;
                    br = ot.bulge_size - bd;
                }

                // Clamp to valid range.
                mm = std::min(mm, max_mm);
                bd = std::min(bd, max_dna);
                br = std::min(br, max_rna);

                if (bd == 0 && br == 0)
                    prof.mm_counts[mm]++;
                else if (bd > 0 && br == 0 && max_dna > 0)
                    prof.dna_counts[mm][bd - 1]++;
                else if (br > 0 && bd == 0 && max_rna > 0)
                    prof.rna_counts[mm][br - 1]++;

                prof.joint_counts[mm][bd][br]++;

                // Per-position mismatch / bulge flags from the aligned strings.
                // target_seq has PAM appended, so first guide_len chars are the guide.
                const std::string& tseq = ot.target_seq;
                const std::string& gseq = ot.guide_seq;

                // Strip PAM from display strings to get the aligned region.
                // pam_at_start → PAM is the first pam_limit chars.
                int g_start = params.pam_at_start ? params.pam_limit : 0;
                int t_start = params.pam_at_start ? params.pam_limit : 0;

                int pos = 0; // guide position counter (skips gap chars)
                for (int i = 0;
                     i < static_cast<int>(tseq.size()) - params.pam_limit && pos < guide_len; ++i)
                {
                    char tc = tseq[t_start + i];
                    char gc =
                        (g_start + i < static_cast<int>(gseq.size())) ? gseq[g_start + i] : 'N';

                    if (tc == '-')
                    {
                        // RNA bulge: gap in target at this guide position.
                        if (pos < guide_len)
                            prof.rna_per_pos[pos]++;
                        pos++;
                    }
                    else if (gc == '-')
                    {
                        // DNA bulge: gap in guide; target has an extra base.
                        if (pos < guide_len)
                            prof.dna_per_pos[pos]++;
                        // do NOT advance pos (this target base has no guide counterpart)
                    }
                    else
                    {
                        // Normal or mismatch position.
                        if (std::islower(static_cast<unsigned char>(tc)))
                        {
                            // Mismatch: record position and nucleotide identity.
                            if (pos < guide_len)
                            {
                                prof.mm_per_pos[pos]++;
                                char uc =
                                    static_cast<char>(std::toupper(static_cast<unsigned char>(tc)));
                                int nuc_idx = -1;
                                if (uc == 'A')
                                    nuc_idx = 0;
                                else if (uc == 'C')
                                    nuc_idx = 1;
                                else if (uc == 'G')
                                    nuc_idx = 2;
                                else if (uc == 'T')
                                    nuc_idx = 3;
                                if (nuc_idx >= 0)
                                    prof.nuc_per_pos[mm][nuc_idx][pos]++;
                            }
                        }
                        pos++;
                    }
                }
            }

            return prof;
        }

        // ---------------------------------------------------------------------------
        // Unique key for deduplication: chr + pos + strand + guide + target.
        // ---------------------------------------------------------------------------
        static std::string hit_key(const OffTarget& ot)
        {
            return ot.chromosome + '\0' + std::to_string(ot.genomic_pos) + '\0' + ot.strand + '\0' +
                   ot.guide_seq + '\0' + ot.target_seq;
        }

    } // anonymous namespace

    // -----------------------------------------------------------------------------
    // load_index
    // -----------------------------------------------------------------------------

    TSTIndex load_index(const std::string& bin_path, int pam_limit, const std::string& chr_name)
    {
        std::ifstream fin(bin_path, std::ios::binary);
        if (!fin.is_open())
            throw std::runtime_error("Cannot open TST index file: " + bin_path);

        TSTIndex idx;
        idx.pam_limit = pam_limit;
        idx.chr_name = chr_name;

        // ---- header ----
        int32_t num_leaves = 0, stored_guide_len = 0;
        fin.read(reinterpret_cast<char*>(&num_leaves), sizeof(int32_t));
        fin.read(reinterpret_cast<char*>(&stored_guide_len), sizeof(int32_t));
        if (!fin)
            throw std::runtime_error("Truncated header in: " + bin_path);

        idx.guide_length = stored_guide_len;
        idx.leaves.resize(num_leaves);

        // ---- leaf array ----
        const int pam_bytes = (pam_limit + 1) / 2;

        for (int i = 0; i < num_leaves; ++i)
        {
            TSTIndexLeaf& leaf = idx.leaves[i];

            int32_t guide_index = 0;
            fin.read(reinterpret_cast<char*>(&guide_index), sizeof(int32_t));
            if (!fin)
                throw std::runtime_error("Truncated leaf in: " + bin_path);
            leaf.guide_index = guide_index;

            leaf.pam_enc.resize(pam_bytes);
            fin.read(reinterpret_cast<char*>(leaf.pam_enc.data()), pam_bytes);
            if (!fin)
                throw std::runtime_error("Truncated PAM in leaf " + std::to_string(i) +
                                         " of: " + bin_path);

            char link_flag = 0;
            fin.get(link_flag);
            if (!fin)
                throw std::runtime_error("Truncated leaf link in: " + bin_path);

            if (link_flag == '0')
            {
                leaf.next = -1;
            }
            else // '_'
            {
                int32_t next_idx = 0;
                fin.read(reinterpret_cast<char*>(&next_idx), sizeof(int32_t));
                if (!fin)
                    throw std::runtime_error("Truncated leaf next in: " + bin_path);
                leaf.next = next_idx;
            }
        }

        // ---- node array ----
        int32_t num_nodes = 0;
        fin.read(reinterpret_cast<char*>(&num_nodes), sizeof(int32_t));
        if (!fin)
            throw std::runtime_error("Truncated node count in: " + bin_path);

        idx.nodes.resize(num_nodes);

        // The TST was serialised in a pre-order recursive manner using the
        // pack_nibbles encoding (two 4-bit IUPAC codes per byte).
        // We deserialise it by reading pairs of nibbles and reconstructing
        // the lo → hi → eq child structure.
        //
        // The sentinel nibble (0xF) signals end-of-sequence ('_') and is
        // followed by a 4-byte leaf index stored verbatim.
        // The null nibble  (0x0) signals "no child" ('0').
        //
        // We replicate the logic of deSerialize() from the legacy code.

        // Read one nibble at a time from the packed byte stream.
        // We buffer at most one unread nibble between calls.
        uint8_t nibble_buf = 0;
        bool nibble_ready = false;

        auto read_nibble = [&]() -> uint8_t
        {
            if (nibble_ready)
            {
                nibble_ready = false;
                return nibble_buf;
            }
            char byte = 0;
            fin.get(byte);
            uint8_t ub = static_cast<uint8_t>(byte);
            nibble_buf = ub & 0x0F;
            nibble_ready = true;
            return (ub >> 4) & 0x0F;
        };

        // Recursive lambda for deserialisation.
        // Returns the number of nodes populated (used only to detect errors).
        std::function<void(int)> deserialise = [&](int cur)
        {
            if (!fin || cur < 0 || cur >= num_nodes)
                return;

            TSTIndexNode& node = idx.nodes[cur];

            // Read splitchar nibble.
            uint8_t sc_nibble = read_nibble();
            node.splitchar_enc = sc_nibble;
            node.splitchar = NUC_DECODE[sc_nibble & 0x0F];

            // lo child
            uint8_t lo_nibble = read_nibble();
            if (lo_nibble == NULL_CHILD_NIBBLE)
            {
                node.lokid = 0;
            }
            else
            {
                // lo_nibble is the splitchar of the lo child's node.
                // Allocate a fresh node index (pre-order: allocate as we go).
                // We need a counter for next available node slot.
                // Because TSTIndex nodes were pre-sized to num_nodes, we use a
                // simple linear scan.  In practice the legacy code fills nodes in
                // pre-order so we just use a running counter passed by ref.
                // We handle this by pushing the nibble back and recursing.
                //
                // Actually: the legacy serialize() writes the splitchar of the
                // child node first, so lo_nibble IS the first nibble of the lo
                // child's record.  We stash it back and recurse.
                //
                // Stash back: next call to read_nibble() must return lo_nibble.
                // We achieve this by noting nibble_ready is false here (we just
                // consumed both nibbles of a byte above).  So we set nibble_buf
                // and nibble_ready directly.
                nibble_buf = lo_nibble;
                nibble_ready = true;

                // Allocate the lo child node slot.
                // We need a mutable counter; capture it by ref via a shared_ptr trick.
                // Simpler: use a static local won't work for re-entrancy.
                // The cleanest solution: pre-allocate all nodes and assign them in
                // DFS order using a shared node counter.
                // We store the counter in a local std::shared_ptr<int> captured by
                // the outer lambda.
                node.lokid = -1; // placeholder; will be fixed below
            }

            // NOTE: The recursive scheme above quickly becomes complex. Instead,
            // we replicate the exact logic from deSerialize() in searchOnTST.cpp
            // which uses a single int pNode passed by reference.  We restart with
            // a cleaner iterative approach below.
            (void)cur; // suppress unused warning in this placeholder path
        };
        (void)deserialise; // suppress unused lambda warning

        // -----------------------------------------------------------------
        // Proper re-implementation of deSerialize() using a node counter.
        // -----------------------------------------------------------------
        // Reset stream to start of node data.
        // (We already read num_nodes above; the stream is positioned right.)
        nibble_ready = false;

        int pNode = 0; // next node index to allocate (root is 0, pre-allocated)

        // Pairs are consumed two nibbles at a time; the legacy code reads one byte
        // = two nibbles, where nibble[0] = first and nibble[1] = second.
        // readPair fills pairNuc_bit[0] and pairNuc_bit[1] from one byte.
        // The flag toggles which nibble of the pair to use next.

        uint8_t pair_byte = 0;
        bool have_pair = false;
        int pair_flag = 1; // 1 = use high nibble next (matches legacy flag=1 at start)

        auto next_nibble_legacy = [&]() -> uint8_t
        {
            if (pair_flag == 1)
            {
                // Read a new byte.
                char b = 0;
                fin.get(b);
                pair_byte = static_cast<uint8_t>(b);
                pair_flag = 0;
                return (pair_byte >> 4) & 0x0F;
            }
            else
            {
                pair_flag = 1;
                return pair_byte & 0x0F;
            }
        };

        // We pre-read the very first nibble before entering deSerialize, just
        // as the legacy code does (it reads the first byte and sets flag=0 before
        // calling deSerialize with flag=0).
        // Read first byte, use high nibble for the root splitchar.
        {
            char b = 0;
            if (fin.get(b))
            {
                pair_byte = static_cast<uint8_t>(b);
                pair_flag = 0; // next call will return low nibble
            }
        }

        std::function<void(int)> deser = [&](int i)
        {
            if (i < 0 || i >= num_nodes || !fin)
                return;

            TSTIndexNode& node = idx.nodes[i];

            // splitchar: already in pair_byte / pair_flag state.
            uint8_t sc;
            if (pair_flag == 0)
            {
                // Use high nibble of current byte.
                sc = (pair_byte >> 4) & 0x0F;
                pair_flag = 1; // next read gives low nibble
            }
            else
            {
                // pair_flag == 1: read new byte, use high nibble.
                char b = 0;
                fin.get(b);
                pair_byte = static_cast<uint8_t>(b);
                sc = (pair_byte >> 4) & 0x0F;
                pair_flag = 0;
            }
            node.splitchar_enc = sc;
            node.splitchar = NUC_DECODE[sc & 0x0F];

            // lo child
            uint8_t lo_nibble = next_nibble_legacy();
            if (lo_nibble != NULL_CHILD_NIBBLE)
            {
                node.lokid = ++pNode;
                // The lo_nibble IS the splitchar of the lo child.
                // We must "push it back" by forcing pair_byte / pair_flag.
                // Simulate: the next call to deser(node.lokid) should consume lo_nibble
                // as its splitchar.  We store lo_nibble in the high-nibble position.
                pair_byte = static_cast<uint8_t>(lo_nibble << 4);
                pair_flag = 0;
                deser(node.lokid);
            }

            // hi child
            uint8_t hi_nibble = next_nibble_legacy();
            if (hi_nibble != NULL_CHILD_NIBBLE)
            {
                node.hikid = ++pNode;
                pair_byte = static_cast<uint8_t>(hi_nibble << 4);
                pair_flag = 0;
                deser(node.hikid);
            }

            // eq child or leaf
            uint8_t eq_nibble = next_nibble_legacy();
            if (eq_nibble == SENTINEL_NIBBLE)
            {
                // Leaf: read the 4-byte leaf index stored verbatim.
                int32_t leaf_idx = 0;
                fin.read(reinterpret_cast<char*>(&leaf_idx), sizeof(int32_t));
                node.eqkid = leaf_idx; // negative value encodes leaf chain
            }
            else
            {
                node.eqkid = ++pNode;
                pair_byte = static_cast<uint8_t>(eq_nibble << 4);
                pair_flag = 0;
                deser(node.eqkid);
            }
        };

        if (num_nodes > 0)
            deser(0);

        fin.close();
        return idx;
    }

    // -----------------------------------------------------------------------------
    // search_offtargets – primary overload
    // -----------------------------------------------------------------------------

    SearchResult search_offtargets(const std::vector<std::string>& bin_paths,
                                   const std::vector<std::string>& guides,
                                   const SearchParams& params)
    {
        if (bin_paths.empty())
            throw std::invalid_argument("bin_paths must not be empty");
        if (guides.empty())
            throw std::invalid_argument("guides must not be empty");
        if (params.max_mismatches < 0)
            throw std::invalid_argument("max_mismatches must be >= 0");
        if (params.max_dna_bulges < 0)
            throw std::invalid_argument("max_dna_bulges must be >= 0");
        if (params.max_rna_bulges < 0)
            throw std::invalid_argument("max_rna_bulges must be >= 0");
        if (params.pam_length <= params.pam_limit)
            throw std::invalid_argument("pam_length must be > pam_limit");
        if (params.pam_limit <= 0)
            throw std::invalid_argument("pam_limit must be > 0");

        const int num_threads = std::max(1, params.num_threads);
        const int n_parts = static_cast<int>(bin_paths.size());
        const int n_guides = static_cast<int>(guides.size());
        const int guide_len = params.guide_length();

        // Pre-encode all guides into 4-bit arrays.
        // Each guide is stored as guide_len uint8_t values.
        std::vector<std::vector<uint8_t>> guide_enc(n_guides, std::vector<uint8_t>(guide_len, 0));
        for (int g = 0; g < n_guides; ++g)
        {
            const std::string& gs = guides[g];
            for (int i = 0; i < guide_len && i < static_cast<int>(gs.size()); ++i)
                guide_enc[g][i] = iupac::encode_genome(gs[i]);
        }

        // Offset for guide length adjustment when the tree was built with a
        // different guide length than what we're searching.  Per legacy code,
        // offset_guide_len = stored_guide_len - search_guide_len.
        // We compute it per-partition after loading.

        // Result accumulator per thread.
        std::vector<std::vector<OffTarget>> thread_hits(num_threads);

        // Per-thread working state.
        std::vector<ThreadState> thread_states(num_threads);
        for (int t = 0; t < num_threads; ++t)
            thread_states[t].init(guide_len, params.max_dna_bulges);

        // Load all partitions up front to allow OpenMP to distribute them.
        std::vector<TSTIndex> indices(n_parts);
        for (int p = 0; p < n_parts; ++p)
        {
            // Derive chr name from the filename: <pam>_<chr>_<part>.bin
            std::string chr_name;
            const std::string& path = bin_paths[p];
            auto slash = path.find_last_of("/\\");
            std::string fname = (slash == std::string::npos) ? path : path.substr(slash + 1);
            // Remove .bin suffix.
            if (fname.size() > 4 && fname.substr(fname.size() - 4) == ".bin")
                fname = fname.substr(0, fname.size() - 4);
            // Format: <pam_seq>_<chr_name>_<part>
            auto first_us = fname.find('_');
            auto last_us = fname.rfind('_');
            if (first_us != std::string::npos && last_us != first_us)
                chr_name = fname.substr(first_us + 1, last_us - first_us - 1);
            else
                chr_name = fname;

            try
            {
                indices[p] = load_index(path, params.pam_limit, chr_name);
            }
            catch (const std::exception& e)
            {
                throw std::runtime_error(std::string("Failed loading index '") + path +
                                         "': " + e.what());
            }
        }

        // Parallel search: distribute (partition × guide) work across threads.
        // We use a coarse-grained schedule: one partition per task, and within
        // each partition we iterate over all guides sequentially.
        // This keeps synchronisation overhead low and avoids false sharing.

#pragma omp parallel num_threads(num_threads)
        {
            int thr = 0;
#ifdef _OPENMP
            thr = omp_get_thread_num();
#endif
            ThreadState& ts = thread_states[thr];

#pragma omp for schedule(dynamic) nowait
            for (int p = 0; p < n_parts; ++p)
            {
                const TSTIndex& idx = indices[p];
                const int offset_glen = idx.guide_length - guide_len;

                for (int g = 0; g < n_guides; ++g)
                {
                    ts.init(guide_len, params.max_dna_bulges);

                    if (idx.nodes.empty())
                        continue;

                    nearsearch(idx,
                               0, // root node
                               guide_enc[g].data(),
                               0, // guide_pos
                               params.max_mismatches, params.max_dna_bulges, params.max_rna_bulges,
                               0, 0, // bulge type counters
                               true, params, idx.chr_name, g, guides[g], offset_glen, ts);

                    // Tag each hit with guide index for profiling.
                    // (guide_idx is stored inside each OffTarget via guide_seq.)
                }

#pragma omp critical
                {
                    thread_hits[thr].insert(thread_hits[thr].end(), ts.hits.begin(), ts.hits.end());
                    ts.hits.clear();
                }
            }
        }

        // Merge all per-thread hit lists.
        std::vector<OffTarget> all_hits;
        for (auto& tv : thread_hits)
            all_hits.insert(all_hits.end(), std::make_move_iterator(tv.begin()),
                            std::make_move_iterator(tv.end()));

        // Deduplicate.
        {
            std::unordered_set<std::string> seen;
            seen.reserve(all_hits.size());
            std::vector<OffTarget> unique_hits;
            unique_hits.reserve(all_hits.size());
            for (auto& ot : all_hits)
            {
                std::string key = hit_key(ot);
                if (seen.insert(std::move(key)).second)
                    unique_hits.push_back(std::move(ot));
            }
            all_hits = std::move(unique_hits);
        }

        // Build per-guide profiles.
        // Group hits by guide_seq (first pam_limit or trailing pam_limit chars
        // are Ns; we match by the guide portion).
        // For simplicity we iterate over all guides and filter hits per guide.
        std::vector<PositionalProfile> profiles;
        profiles.reserve(n_guides);

        for (int g = 0; g < n_guides; ++g)
        {
            // Build PAM-appended guide string for matching.
            std::string g_full;
            if (params.pam_at_start)
                g_full = std::string(params.pam_limit, 'N') + guides[g];
            else
                g_full = guides[g] + std::string(params.pam_limit, 'N');

            std::vector<OffTarget> guide_hits;
            for (const OffTarget& ot : all_hits)
                if (ot.guide_seq == g_full)
                    guide_hits.push_back(ot);

            profiles.push_back(build_profile(guides[g], guide_hits, params));
        }

        SearchResult result;
        result.off_targets = std::move(all_hits);
        result.profiles = std::move(profiles);
        return result;
    }

    // -----------------------------------------------------------------------------
    // search_offtargets – scalar convenience overload
    // -----------------------------------------------------------------------------

    SearchResult search_offtargets(const std::vector<std::string>& bin_paths,
                                   const std::vector<std::string>& guides, int max_mismatches,
                                   int max_dna_bulges, int max_rna_bulges, int pam_length,
                                   int pam_limit, bool pam_at_start, int num_threads)
    {
        SearchParams p;
        p.max_mismatches = max_mismatches;
        p.max_dna_bulges = max_dna_bulges;
        p.max_rna_bulges = max_rna_bulges;
        p.pam_length = pam_length;
        p.pam_limit = pam_limit;
        p.pam_at_start = pam_at_start;
        p.num_threads = std::max(1, num_threads);
        return search_offtargets(bin_paths, guides, p);
    }

} // namespace crispritz