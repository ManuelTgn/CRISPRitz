#pragma once

#include "offtarget.hpp"
#include "tst_utils.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace crispritz
{

    // =============================================================================
    // SearchParams
    // =============================================================================

    /**
     * @brief Runtime parameters controlling one off-target search run.
     *
     * All limits are inclusive: max_mismatches = 2 allows 0, 1, or 2 substitutions.
     * DNA bulge  = gap in the guide  / extra base in the genomic target.
     * RNA bulge  = gap in the target / extra base in the guide.
     *
     * @param max_mismatches  Maximum number of substitution mismatches allowed.
     * @param max_dna_bulges  Maximum number of DNA bulges (insertions in target).
     * @param max_rna_bulges  Maximum number of RNA bulges (deletions in target).
     * @param pam_length      Total length of the PAM + guide pattern
     *                        (e.g. 23 for NNNNNNNNNNNNNNNNNNNNNGG).
     * @param pam_limit       Length of the PAM-only portion (e.g. 3 for NGG).
     * @param pam_at_start    True when PAM precedes the guide (e.g. Cas12a TTTN).
     * @param num_threads     Number of OpenMP threads for the parallel search.
     */
    struct SearchParams
    {
        int max_mismatches = 0;
        int max_dna_bulges = 0;
        int max_rna_bulges = 0;
        int pam_length = 0;
        int pam_limit = 0;
        bool pam_at_start = false;
        int num_threads = 1;

        /** @return Guide length derived as pam_length - pam_limit. */
        int guide_length() const noexcept
        {
            return pam_length - pam_limit;
        }
    };

    // =============================================================================
    // PositionalProfile
    // =============================================================================

    /**
     * @brief Per-guide positional mismatch and bulge counts accumulated across
     *        all searched chromosomes.
     *
     * Mirrors the legacy profiling matrices written to `.profile.xls`,
     * `.profile_dna.xls`, `.profile_rna.xls`, and `.profile_complete.xls`.
     *
     * All index-by-position vectors have length == guide_len.
     * Mismatch-count histograms have length == max_mm + 1.
     */
    struct PositionalProfile
    {
        std::string guide; ///< Guide sequence (with PAM Ns appended).

        // -- per-position counts ------------------------------------------------

        /** @brief mm_per_pos[i]  = total mismatches seen at guide position i. */
        std::vector<int> mm_per_pos;

        /** @brief dna_per_pos[i] = total DNA-bulge events at guide position i. */
        std::vector<int> dna_per_pos;

        /** @brief rna_per_pos[i] = total RNA-bulge events at guide position i. */
        std::vector<int> rna_per_pos;

        // -- summary histograms -------------------------------------------------

        /**
         * @brief mm_counts[k] = number of targets with exactly k mismatches
         *        (and 0 bulges).  Length = max_mm + 1.
         */
        std::vector<int> mm_counts;

        /**
         * @brief dna_counts[k][j] = targets with k mismatches and j+1 DNA bulges.
         *        Dimensions: (max_mm+1) x max_dna.
         */
        std::vector<std::vector<int>> dna_counts;

        /**
         * @brief rna_counts[k][j] = targets with k mismatches and j+1 RNA bulges.
         *        Dimensions: (max_mm+1) x max_rna.
         */
        std::vector<std::vector<int>> rna_counts;

        /**
         * @brief joint_counts[k][bd][br] = targets with k mismatches, bd DNA
         *        bulges, and br RNA bulges.
         *        Dimensions: (max_mm+1) x (max_dna+1) x (max_rna+1).
         */
        std::vector<std::vector<std::vector<int>>> joint_counts;

        /**
         * @brief nuc_per_pos[k][n][i] = count of nucleotide n at position i in
         *        targets with k mismatches (n: 0=A, 1=C, 2=G, 3=T).
         *        Dimensions: (max_mm+1) x 4 x guide_len.
         */
        std::vector<std::vector<std::vector<int>>> nuc_per_pos;
    };

    // =============================================================================
    // SearchResult
    // =============================================================================

    /**
     * @brief Aggregated output of a complete off-target search run.
     *
     * @param off_targets  All individual hits, deduplicated, one per row of the
     *                     legacy `.targets.txt` output.
     * @param profiles     One PositionalProfile per input guide, in input order.
     */
    struct SearchResult
    {
        std::vector<OffTarget> off_targets;
        std::vector<PositionalProfile> profiles;
    };

    // =============================================================================
    // In-memory TST index (loaded from .bin files)
    // =============================================================================

    /**
     * @brief A single node in a loaded TST partition.
     *
     * Matches the legacy `Tnode` struct from searchOnTST.cpp.
     * eqkid < 0 encodes a leaf: leaf_index = (-eqkid - 1).
     */
    struct TSTIndexNode
    {
        char splitchar = '\0';     ///< IUPAC character stored at this node.
        uint8_t splitchar_enc = 0; ///< 4-bit IUPAC encoding (for fast AND match).
        int lokid = 0;             ///< Left child index.
        int hikid = 0;             ///< Right child index.
        int eqkid = 0;             ///< Equal child, or -(leaf_idx+1) when leaf.
    };

    /**
     * @brief Metadata for one PAM site stored at a TST leaf.
     *
     * Multiple genomic positions sharing the same guide sequence are chained
     * via the @p next field (linked list, -1 = end of chain).
     *
     * @param guide_index  Signed genomic position: positive = forward strand,
     *                     negative = reverse strand (take absolute value).
     * @param pam_enc      Bit-packed PAM bytes (two IUPAC nibbles per byte,
     *                     high nibble first), length = ceil(pam_limit / 2).
     * @param next         Index of next leaf with the same guide_seq, or -1.
     */
    struct TSTIndexLeaf
    {
        int guide_index = 0;
        std::vector<uint8_t> pam_enc;
        int next = -1;
    };

    /**
     * @brief A fully loaded TST partition ready for off-target search.
     *
     * @param nodes        All TST internal nodes (root is index 0).
     * @param leaves       All leaf entries (genomic positions + packed PAM).
     * @param guide_length Guide length stored in the .bin header.
     * @param pam_limit    PAM-only length (bytes per leaf = ceil(pam_limit/2)).
     * @param chr_name     Chromosome / contig this partition belongs to.
     */
    struct TSTIndex
    {
        std::vector<TSTIndexNode> nodes;
        std::vector<TSTIndexLeaf> leaves;
        int guide_length = 0;
        int pam_limit = 0;
        std::string chr_name;
    };

    // =============================================================================
    // Public free functions
    // =============================================================================

    /**
     * @brief Load a single `.bin` TST partition file into a TSTIndex.
     *
     * Binary layout (byte-exact with TernarySearchTree::write_partition):
     * @code
     *   [int32]               num_leaves
     *   [int32]               guide_length
     *   for each leaf:
     *     [int32]             guide_index  (signed genomic pos)
     *     [ceil(pam_limit/2)] packed PAM nibbles (high nibble first)
     *     ['0']               no next leaf, OR
     *     ['_' + int32]       next leaf index
     *   [int32]               num_nodes
     *   [variable]            serialised TST nodes (pack_nibbles encoding)
     * @endcode
     *
     * @param bin_path   Path to the `.bin` partition file.
     * @param pam_limit  PAM-only length (determines bytes read per leaf).
     * @param chr_name   Chromosome name stored in the returned index.
     * @return           A fully populated TSTIndex.
     *
     * @throws std::runtime_error on any I/O or format error.
     */
    TSTIndex load_index(const std::string& bin_path, int pam_limit, const std::string& chr_name);

    /**
     * @brief Search all provided TST partitions for off-target sites of every
     *        input guide, returning hits and per-guide positional profiles.
     *
     * Algorithm overview:
     *  1. Load each `.bin` file into a TSTIndex (parallelised with OpenMP).
     *  2. For each (chromosome-partition, guide) pair run the recursive TST
     *     near-search tolerating up to max_mismatches substitutions,
     *     max_dna_bulges gaps in the guide, and max_rna_bulges gaps in the target.
     *  3. Merge per-thread hit lists, deduplicate identical target strings
     *     (same guide_seq + target_seq + chromosome + position).
     *  4. Accumulate per-guide PositionalProfile from the deduplicated hits.
     *
     * @param bin_paths  Paths to `.bin` partition files (all chromosomes).
     * @param guides     Guide sequences WITHOUT PAM, one per query.
     * @param params     Search configuration (limits, PAM info, threading).
     * @return           SearchResult with all hits and positional profiles.
     *
     * @throws std::invalid_argument if any limit is negative or guides is empty.
     * @throws std::runtime_error    on I/O failure or malformed .bin file.
     *
     * @complexity O(G × P × L^(M+B)) where G = #guides, P = #partitions,
     *             L = guide length, M = max_mismatches, B = max total bulges.
     */
    SearchResult search_offtargets(const std::vector<std::string>& bin_paths,
                                   const std::vector<std::string>& guides,
                                   const SearchParams& params);

    /**
     * @brief Scalar-parameter convenience overload; exposed directly to Python.
     *
     * Constructs a SearchParams internally and delegates to the primary overload.
     *
     * @param bin_paths       Paths to `.bin` partition files.
     * @param guides          Guide sequences WITHOUT PAM.
     * @param max_mismatches  Maximum substitutions allowed (>= 0).
     * @param max_dna_bulges  Maximum DNA bulges allowed  (>= 0).
     * @param max_rna_bulges  Maximum RNA bulges allowed  (>= 0).
     * @param pam_length      Full PAM+guide pattern length.
     * @param pam_limit       PAM-only length.
     * @param pam_at_start    True for PAM-upstream systems (e.g. Cas12a).
     * @param num_threads     OpenMP thread count (clamped to 1 if <= 0).
     * @return                SearchResult with hits and positional profiles.
     *
     * @throws std::invalid_argument if any limit is negative.
     * @throws std::runtime_error    on I/O or search failure.
     */
    SearchResult search_offtargets(const std::vector<std::string>& bin_paths,
                                   const std::vector<std::string>& guides, int max_mismatches,
                                   int max_dna_bulges, int max_rna_bulges, int pam_length,
                                   int pam_limit, bool pam_at_start, int num_threads = 1);

} // namespace crispritz
