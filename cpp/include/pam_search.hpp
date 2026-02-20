#ifndef PAM_SEARCH_HPP
#define PAM_SEARCH_HPP

#include <string>
#include <string_view>
#include <vector>
#include <cstdint>

namespace pam {

struct SearchParams {
    int pam_length;
    int pam_limit;           // How much of the PAM to check
    int max_mismatches;
    int max_bulges;
    bool pam_at_start;       // PAM at 3' (start) vs 5' (end)
    
    SearchParams(int len, int limit, int mm, int bulges, bool at_start)
        : pam_length(len), pam_limit(limit), max_mismatches(mm), 
          max_bulges(bulges), pam_at_start(at_start) {}
};

// Compact genome representation using 2 bits per nucleotide packed in bytes
class CompactGenome {
public:
    explicit CompactGenome(std::string_view sequence);
    
    uint8_t operator[](size_t pos) const {
        return (data_[pos >> 1] >> ((pos & 1) << 2)) & 0x0F;
    }
    
    size_t size() const { return length_; }
    
    // Allow access to raw data for optimized searching
    const uint8_t* data() const noexcept { return data_.data(); }

    size_t bytes() const noexcept { return data_.size(); }
    
    // Friend declaration for search functions that need direct access
    friend std::vector<int> search_pam_sites_fast(
        std::string_view pam_sequence,
        const CompactGenome& genome,
        const SearchParams& params
    );
    
private:
    std::vector<uint8_t> data_;  // Packed 4-bit encodings (2 per byte)
    size_t length_;
};

// Main PAM search function
// Returns positions where PAM matches (positive for forward, negative for reverse)
std::vector<int> search_pam_sites(
    std::string_view pam_sequence,
    std::string_view genome_sequence,
    const SearchParams& params
);

// Optimized version using pre-encoded genome
std::vector<int> search_pam_sites_fast(
    std::string_view pam_sequence,
    const CompactGenome& genome,
    const SearchParams& params
);

} // namespace pam

#endif // PAM_SEARCH_HPP