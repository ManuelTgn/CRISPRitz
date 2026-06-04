#pragma once

#include <string>
#include <stdexcept>

namespace crispritz {

/**
 * @brief Immutable configuration parameters for an off-target search execution.
 * * Bundling these parameters guarantees that invariant search constraints 
 * remain unmodified across recursive backtracking branches.
 */
struct SearchConfig {
    int max_mismatches = 0;
    int max_dna_bulges = 0;
    int max_rna_bulges = 0;
    
    int guide_length = 0;
    bool pam_at_start = false;
    std::string pam_sequence;

    /**
     * @brief Validates the configuration boundaries before starting a search.
     * @throws std::invalid_argument if parameters are logically inconsistent.
     */
    void validate() const {
        if (max_mismatches < 0 || max_dna_bulges < 0 || max_rna_bulges < 0) {
            throw std::invalid_argument("Search limits (mismatches, bulges) cannot be negative.");
        }
        if (guide_length <= 0) {
            throw std::invalid_argument("Guide length must be strictly positive.");
        }
        if (pam_sequence.empty()) {
            throw std::invalid_argument("PAM sequence cannot be empty.");
        }
    }
};

} // namespace crispritz