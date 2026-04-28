#pragma once

/**
 * @file offtargets.hpp
 * @brief Base off-targets representation in crispritz.
 * 
 * Defines the per-hit result record used to accumulate off-target matches 
 * before writing output.
 */

#include <string>

namespace crispritz {

// -------------------------------------------------------------------------
// Per-hit result record
// -------------------------------------------------------------------------

/**
 * @brief Holds all output fields for a single off-target hit.
 *
 * Instances are accumulated per thread and flushed to the output file under
 * a critical section.  Using a value struct avoids pointer chasing and makes
 * copying trivially correct.
 */
struct OffTargetHit
{
    std::string guide_str;        ///< Guide sequence with PAM placeholders
    std::string target_str;       ///< Target sequence (mismatches lowercase)
    std::string chrom;            ///< Chromosome name
    int position{0};      ///< Genomic position (1-based, strand-adjusted)
    int cluster_pos{0};   ///< Position ignoring bulge offset
    char direction{'+'};   ///< Strand direction ('+' or '-')
    int mismatches{0};    ///< Number of mismatches
    int bulge_size{0};    ///< Total bulge length
    std::string bulge_type{"X"};  ///< "X" | "DNA" | "RNA" | "RNA,DNA"
};

}  // namespace crispritz

