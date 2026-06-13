#include "nucleotide_encoding.hpp"

#include <algorithm>

namespace pam {
std::string reverse_complement(std::string_view seq) {
  std::string result;
  result.reserve(seq.length());

  // Iterate in reverse and complement each base
  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    result += NucleotideEncoder::complement(*it);
  }

  return result;
}

} // namespace pam
