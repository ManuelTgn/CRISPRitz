#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

namespace crispritz
{

    // =========================================================================
    // TNode
    // =========================================================================

    /**
     * @brief A single internal node of the Ternary Search Tree.
     *
     * Indices are used instead of raw pointers so that the node array can be
     * stored in a contiguous std::vector without invalidating references on
     * reallocation.  A value of 0 in any child field means "no child" (the root
     * is always at index 0, so child index 0 is reserved as the null sentinel).
     *
     * The eqkid field doubles as a leaf pointer: when negative, the node IS a
     * leaf and the leaf index is @c (-eqkid - 1), matching the legacy encoding
     * used by the binary .bin format and searchOnTST.cpp.
     *
     * @note splitchar_enc stores the 4-bit IUPAC genome encoding
     *       (iupac::encode_genome) of the node character so that match checks
     *       during search reduce to a single bitwise AND, avoiding a second
     *       table lookup.
     * 
     * @param splitchar     IUPAC character stored in the current node. Ambiguous
     *                      nucleotides are denoted by the extended alphabet 
     *                      characters.
     * 
     * @param lokid         Number of the left children of the current node. This
     *                      value denotes how many children store a character
     *                      that is lexicographically lower than the character 
     *                      stored in the current node.
     * 
     * @param hikid         Number of the right children of the current node. This
     *                      value denotes how many children store a character
     *                      that is lexicographically higher than the character 
     *                      stored in the current node.
     * 
     * @param eqkid         Number of the equal children of the current node. This
     *                      value denotes how many children store a character
     *                      that is lexicographically equal to the character 
     *                      stored in the current node.
     */
    struct TNode
    {
        char    splitchar;           ///< IUPAC character stored at this node.
        int     lokid         = 0;   ///< Index of the left  (less-than)    child node.
        int     hikid         = 0;   ///< Index of the right (greater-than) child node.
        int     eqkid         = 0;   ///< Index of the equal child, or -(leaf_idx+1).
    };

    // =========================================================================
    // TLeaf
    // =========================================================================

    /**
     * @brief Metadata stored at a TST leaf, one per PAM site extracted from the
     *        genome.
     *
     * The legacy code stored leaves in a parallel array indexed via the negative
     * eqkid values of terminal TSTNodes.  That layout is preserved here so the
     * binary serialization format does not change.
     *
     * @note The two legacy fields @c pam_dna and @c pam_dna_bit that appeared in
     *       earlier revisions have been removed; they were never written to disk
     *       and inflated struct size by ~40 bytes per leaf (tens of MiB on a
     *       whole-genome index).  All PAM data is carried exclusively in
     *       @c pam_seq_enc.
     *
     * @param guide_index  Signed genomic position.  Positive values denote the
     *                     forward strand; negative values denote the reverse
     *                     strand (convention shared with searchOnTST.cpp).
     * @param guide_seq    The guide-length subsequence extracted from the genome,
     *                     stored 5'→3' AFTER any reversal/RC so that lexicographic
     *                     sort and TST insertion operate on the canonical form used
     *                     during search.  Contains only plain characters A/C/G/T —
     *                     never encoded nibble values.
     * @param pam_seq_enc  Bit-packed PAM nucleotides (two @c encode_pam nibbles per
     *                     byte, high nibble first) exactly as written to the .bin
     *                     leaf section.
     * @param next         Index of the next leaf sharing the same guide_seq (linked
     *                     list of genomic positions for the same k-mer).  0 means
     *                     no next leaf.
     */
    struct TLeaf
    {
        int                  guide_index;       ///< Signed genomic position.
        std::string          guide_dna;         ///< Plain-char guide subsequence.
        std::string          pam_dna;           ///< Packed PAM bytes (encode_pam).
        int                  next   = 0;        ///< Next leaf with same guide_seq.
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
     *   TernarySearchTree tst(sequence, chr_name, pam_seq, pam_length,
     *                         pam_limit, pam_at_start, outdir,
     *                         max_bulges, num_threads);
     *   tst.build();   // PAM search → extraction → sort → insert → save
     * @endcode
     *
     * The serialized .bin files are consumed by the search module without
     * modification; the format is byte-for-byte identical to the legacy
     * mainParallel.cpp output.
     *
     * Thread safety: build() uses OpenMP to parallelize the PAM search phase.
     * All other methods are single-threaded and must not be called concurrently.
     */
    class TernarySearchTree
    {
      public:

        /**
         * @brief Construct the TST builder.
         *
         * No computation is performed here; call build() to run the pipeline.
         *
         * @param sequence     Full genomic sequence (single chromosome, uppercase
         *                     IUPAC characters).
         * @param chr_name     Chromosome / contig identifier used in output filenames.
         * @param pam_seq      The PAM-only string (e.g. @c "NGG"), without guide Ns.
         * @param pam_length   Total length of the PAM+guide pattern (e.g. 23 for
         *                     @c NNNNNNNNNNNNNNNNNNNNNGG).
         * @param pam_limit    Length of the PAM portion only (e.g. 3 for @c "NGG").
         * @param pam_at_start True when the PAM precedes the guide (e.g. Cas12a).
         * @param max_bulges   Maximum number of bulges; extra bases are extracted per
         *                     site to allow bulge-aware search later.
         * @param num_threads  Number of OpenMP threads for the PAM search phase.
         *
         * @throws std::runtime_error if guide_length ≤ 0 or pam_limit ≤ 0 or
         *         outdir is empty.
         */
        TernarySearchTree(const std::string& sequence, const std::string& chr_name,
                          const std::string& pam_seq, int pam_length, int pam_limit,
                          bool pam_at_start, int max_bulges = 0, int num_threads = 1);

        /**
         * @brief Run the full build pipeline:
         *        1. PAM search  → site positions
         *        2. Sequence extraction → TSTLeaf array (N-containing windows dropped)
         *        3. Lexicographic sort on guide_seq
         *        4. Balanced median insertion into TST nodes
         *        5. Serialization to .bin partition file(s) under outdir_
         *
         * @throws std::runtime_error if all PAM sites are discarded (N-filtering),
         *         or if a .bin file cannot be opened for writing.
         */
        void build();

        /**
         * @brief Write the constructed TST to disk.
         *
         * Large chromosomes are split into chunks of LEAVES_PER_GROUP leaves,
         * each written to a separate .bin file named:
         * @code
         *   <outdir_>/<pam_seq_>_<chr_name_>_<part>.bin
         * @endcode
         *
         * This method is called internally by build() and is exposed for testing.
         * 
         * @param outdir    The target output directory where binary partition 
         *                  files will be stored.
         *
         * @throws std::runtime_error on any I/O failure.
         */
        void save(const std::string& outdir);

      private:

        // ------------------------------------------------------------------
        // Configuration (set in constructor, immutable after that)
        // ------------------------------------------------------------------
        std::string sequence_;
        std::string chr_name_;
        std::string pam_seq_;
        int         pam_length_;
        int         pam_limit_;
        int         guide_length_; ///< Derived: pam_length_ - pam_limit_.
        bool        pam_at_start_;
        std::string outdir_;
        int         max_bulges_;
        int         num_threads_;
        std::string pam_rna_;

        /** @brief Maximum number of leaves stored in a single .bin partition. */
        static constexpr int LEAVES_PER_GROUP = 5'000'000;

        // ---------------------------------------------------------------------
        // Build state (mutated by build() and save())
        // ---------------------------------------------------------------------
        std::vector<TLeaf> leaves_;   ///< Extracted and sorted target leaves.

        // ---------------------------------------------------------------------
        // TST insertion helpers
        // ---------------------------------------------------------------------
        void insert_target(const std::string& target, int global_dx, int relative_idx,
                            std::vector<TNode>& tree, int& node_used);
        void build_sub_tree(int l, int r, std::vector<TNode>& tree, int& node_used,
                            int start);

        // ---------------------------------------------------------------------
        // Serialization helpers
        // ---------------------------------------------------------------------
        void write_pair_nuc(std::ostream& os, char pair_nuc[2], uint8_t& bit_nuc) const;
        bool update_pair_nuc(char pair_nuc[2], bool switch_node, char nt, std::ostream& os, uint8_t& bit_nuc) const;
        void serialize_node(int node_idx, const std::vector<TNode>& tree, std::ostream& os, 
                            bool& switch_node, char pair_nuc[2], uint8_t& bit_nuc) const;
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
     * @param sequence     Full genomic sequence (single chromosome, uppercase).
     * @param chr_name     Chromosome name used in output filenames.
     * @param pam_seq      PAM-only string (e.g. @c "NGG").
     * @param pam_length   Full pattern length including guide Ns.
     * @param pam_limit    PAM-only length.
     * @param pam_at_start True for PAM-upstream systems (e.g. Cas12a).
     * @param outdir       Directory where .bin files will be written (must exist).
     * @param max_bulges   Maximum bulge count (default 0).
     * @param num_threads  Number of OpenMP threads for PAM search (default 1).
     *
     * @throws std::runtime_error on any build or I/O failure.
     */
    inline void build_tree(const std::string& sequence, const std::string& chr_name,
                    const std::string& pam_seq, int pam_length, int pam_limit,
                    bool pam_at_start, const std::string& outdir,
                    int max_bulges = 0, int num_threads = 1)
    {
        try
        {
            TernarySearchTree tst(sequence, chr_name, pam_seq, pam_length, pam_limit, 
                                    pam_at_start, max_bulges, num_threads);
            tst.build();
            tst.save(outdir);
        }
        catch(const std::exception& e)
        {
            throw std::runtime_error(std::string("Failed building Ternary Search Tree: ") + e.what());
        }
    }

} // namespace crispritz