#pragma once

#include <string>
#include <cstdint>

namespace crispritz {

/**
 * @brief Represents the genomic strand orientation.
 */
enum class Strand : char {
    Positive = '+',
    Negative = '-'
};

/**
 * @brief Categorizes the type of structural variation in a match.
 */
enum class BulgeType : uint8_t {
    None,
    DNA,  // Extra bases in the target genome (Insertion relative to guide)
    RNA,  // Extra bases in the guide RNA (Deletion relative to genome)
    Both  // Complex matches containing both DNA and RNA bulges
};

/**
 * @brief Encapsulates a resolved off-target match against the genome.
 * * This structure holds all structural and positional data necessary for 
 * downstream reporting and profiling matrices. It is completely decoupled 
 * from the TST traversal mechanisms.
 */
struct OffTargetMatch {
    // Owning strings: necessary because traversal buffers are transient stack arrays.
    std::string guide_sequence;       
    std::string target_sequence;      
    
    int mismatches = 0;               
    int dna_bulges = 0;               
    int rna_bulges = 0;               
    
    BulgeType bulge_type = BulgeType::None; 
    
    uint64_t genomic_position = 0;    // Absolute position on the chromosome
    uint64_t cluster_position = 0;    // Cluster coordinate (legacy reporting compatibility)
    Strand strand_direction = Strand::Positive;

    /**
     * @brief Determines the total combined bulge size.
     */
    [[nodiscard]] int total_bulges() const noexcept {
        return dna_bulges + rna_bulges;
    }
    
    /**
     * @brief Helper to resolve the structural bulge classification.
     */
    [[nodiscard]] static constexpr BulgeType compute_bulge_type(int dna_bul, int rna_bul) noexcept {
        if (dna_bul > 0 && rna_bul > 0) return BulgeType::Both;
        if (dna_bul > 0) return BulgeType::DNA;
        if (rna_bul > 0) return BulgeType::RNA;
        return BulgeType::None;
    }
};

} // namespace crispritz