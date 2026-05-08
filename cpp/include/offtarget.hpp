#ifndef OFFTARGET_HPP
#define OFFTARGET_HPP

#include <string>

namespace crispritz
{

    // =============================================================================
    // OffTarget
    // =============================================================================

    /**
     * @brief A single off-target hit produced by the TST near-search.
     *
     * Column layout mirrors the legacy `.targets.txt` output so downstream
     * Python code can consume it without adaptation:
     *   bulge_type  crRNA  DNA  chromosome  position  cluster_pos
     *   strand  mismatches  bulge_size  total
     *
     * @note genomic_pos is 1-based. The sign is resolved to the strand field
     *       before storage so genomic_pos is always positive here.
     */
    struct OffTarget
    {
        std::string bulge_type; ///< "X", "DNA", "RNA", or "DNA,RNA".
        std::string guide_seq;  ///< Guide sequence with PAM N-placeholder appended.
        std::string target_seq; ///< Genomic target (mismatch bases in lower-case).
        std::string chromosome; ///< Contig / chromosome name.
        int genomic_pos = 0;    ///< 1-based genomic start position.
        int cluster_pos = 0;    ///< Cluster-normalised start position.
        char strand = '+';      ///< '+' (forward) or '-' (reverse).
        int mismatches = 0;     ///< Number of substitutions.
        int bulge_size = 0;     ///< Total bulge bases (DNA + RNA).
        int total_score = 0;    ///< mismatches + bulge_size.
    };

} // namespace crispritz

#endif // OFFTARGET_HPP