#ifndef PAM_SEARCH_HPP
#define PAM_SEARCH_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace pam {

struct SearchParams {
  int pam_length;
  int pam_limit;     // How much of the PAM to check
  bool pam_at_start; // PAM at 3' (start) vs 5' (end)
  int num_threads;

  SearchParams(int len, int limit, bool at_start, int threads = 1)
      : pam_length(len), pam_limit(limit), pam_at_start(at_start),
        num_threads(threads > 0 ? threads : 1) {}
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
  const uint8_t *data() const noexcept { return data_.data(); }

  size_t bytes() const noexcept { return data_.size(); }

  // Friend declaration for search functions that need direct access
  friend std::vector<int> search_pam_sites_fast(std::string_view pam_sequence,
                                                const CompactGenome &genome,
                                                const SearchParams &params);

private:
  std::vector<uint8_t> data_; // Packed 4-bit encodings (2 per byte)
  size_t length_;
};

// Main PAM search function
// Returns positions where PAM matches (positive for forward, negative for
// reverse)
std::vector<int> search_pam_sites(std::string_view pam_sequence,
                                  std::string_view genome_sequence,
                                  const SearchParams &params);

// Optimized version using pre-encoded genome
std::vector<int> search_pam_sites_fast(std::string_view pam_sequence,
                                       const CompactGenome &genome,
                                       const SearchParams &params);

} // namespace pam

#endif // PAM_SEARCH_HPP