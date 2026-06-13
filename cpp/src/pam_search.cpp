#include "pam_search.hpp"

#include "nucleotide_encoding.hpp"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace pam {

CompactGenome::CompactGenome(std::string_view sequence)
    : length_(sequence.length()) {
  // Allocate space for packed data (2 nucleotides per byte)
  data_.resize((length_ + 1) / 2);

  for (size_t i = 0; i < length_; ++i) {
    uint8_t encoded = NucleotideEncoder::encode_genome(sequence[i]);
    size_t byte_idx = i >> 1;
    size_t shift = (i & 1) << 2; // 0 or 4
    data_[byte_idx] |= (encoded << shift);
  }
}

namespace {

// Check if PAM matches at given position with mismatch tolerance
inline bool check_pam_match(const uint8_t *genome_data, const uint8_t *pam_data,
                            size_t pos, int pam_limit) {

  for (int i = 0; i < pam_limit; ++i) {
    uint8_t g_byte = genome_data[(pos + i) >> 1];
    uint8_t g_shift = ((pos + i) & 1) << 2;
    uint8_t genome_enc = (g_byte >> g_shift) & 0x0F;

    uint8_t p_byte = pam_data[i >> 1];
    uint8_t p_shift = (i & 1) << 2;
    uint8_t pam_enc = (p_byte >> p_shift) & 0x0F;

    // Bitwise AND checks if nucleotides match (IUPAC ambiguity codes)
    if ((genome_enc & pam_enc) == 0) {
      return false;
    }
  }

  return true;
}

// Encode PAM sequence into compact format
std::vector<uint8_t> encode_pam(std::string_view pam, bool reverse_comp) {
  std::string pam_str{pam};
  if (reverse_comp) {
    pam_str = reverse_complement(pam);
  }

  std::vector<uint8_t> encoded((pam_str.length() + 1) / 2, 0);

  for (size_t i = 0; i < pam_str.length(); ++i) {
    uint8_t enc = NucleotideEncoder::encode_pam(pam_str[i]);
    size_t byte_idx = i >> 1;
    size_t shift = (i & 1) << 2;
    encoded[byte_idx] |= (enc << shift);
  }

  return encoded;
}

} // anonymous namespace

std::vector<int> search_pam_sites_fast(std::string_view pam_sequence,
                                       const CompactGenome &genome,
                                       const SearchParams &params) {
  std::vector<int> indices;
  indices.reserve(10000); // Reserve space to reduce reallocations

  const size_t genome_len = genome.size();
  const int pam_len = params.pam_length;
  const int pam_limit = params.pam_limit;
  const bool pam_at_start = params.pam_at_start;

  // Pre-encode PAM sequences
  auto pam_fwd = encode_pam(pam_sequence, false);
  auto pam_rev = encode_pam(pam_sequence, true);

  const uint8_t *pam_fwd_data = pam_fwd.data();
  const uint8_t *pam_rev_data = pam_rev.data();

  // Calculate search bounds
  const size_t search_end = genome_len - pam_limit;

  if (!pam_at_start) { // PAM at 5' (end of guide)
#pragma omp parallel num_threads(params.num_threads)
    {
      std::vector<int> local_indices;
      local_indices.reserve(1000);

#pragma omp for schedule(static) nowait
      for (size_t pos = 0; pos < search_end; ++pos) {
        // Check forward strand
        if (check_pam_match(genome.data(), pam_fwd_data, pos, pam_limit)) {
          int guide_start =
              static_cast<int>(pos + pam_limit - 1) - (pam_len - 1);
          if (guide_start >= 0) {
            local_indices.push_back(guide_start);
          }
        }

        // Check reverse strand
        if (check_pam_match(genome.data(), pam_rev_data, pos, pam_limit)) {
          if (pos <= genome_len - pam_len) {
            local_indices.push_back(-static_cast<int>(pos));
          }
        }
      }

#pragma omp critical
      {
        indices.insert(indices.end(), local_indices.begin(),
                       local_indices.end());
      }
    }
  } else { // PAM at 3' (start of guide)
#pragma omp parallel num_threads(params.num_threads)
    {
      std::vector<int> local_indices;
      local_indices.reserve(1000);

#pragma omp for schedule(static) nowait
      for (size_t pos = 0; pos < search_end; ++pos) {
        // Check forward strand
        if (check_pam_match(genome.data(), pam_fwd_data, pos, pam_limit)) {
          if (pos <= genome_len - pam_len) {
            local_indices.push_back(-static_cast<int>(pos));
          }
        }

        // Check reverse strand
        if (check_pam_match(genome.data(), pam_rev_data, pos, pam_limit)) {
          int guide_start =
              static_cast<int>(pos + pam_limit - 1) - (pam_len - 1);
          if (guide_start >= 0) {
            local_indices.push_back(guide_start);
          }
        }
      }

#pragma omp critical
      {
        indices.insert(indices.end(), local_indices.begin(),
                       local_indices.end());
      }
    }
  }

  return indices;
}

std::vector<int> search_pam_sites(std::string_view pam_sequence,
                                  std::string_view genome_sequence,
                                  const SearchParams &params) {
  CompactGenome genome(genome_sequence);
  return search_pam_sites_fast(pam_sequence, genome, params);
}

} // namespace pam