#ifndef NUCLEOTIDE_ENCODING_HPP
#define NUCLEOTIDE_ENCODING_HPP

#include <cstdint>
#include <string>
#include <string_view>

namespace pam {

// Efficient 4-bit encoding for IUPAC nucleotides
// A=0001, C=0010, G=0100, T=1000, N=1111 (all), etc.
class NucleotideEncoder {
public:
  // Encoding tables as constexpr for compile-time optimization
  static constexpr uint8_t encode_genome(char c) {
    switch (c) {
    case 'A':
      return 0b0001;
    case 'C':
      return 0b0010;
    case 'G':
      return 0b0100;
    case 'T':
      return 0b1000;
    case 'R':
      return 0b0101; // A or G
    case 'Y':
      return 0b1010; // C or T
    case 'S':
      return 0b0110; // G or C
    case 'W':
      return 0b1001; // A or T
    case 'K':
      return 0b1100; // G or T
    case 'M':
      return 0b0011; // A or C
    case 'B':
      return 0b1110; // C or G or T
    case 'D':
      return 0b1101; // A or G or T
    case 'H':
      return 0b1011; // A or C or T
    case 'V':
      return 0b0111; // A or C or G
    default:
      return 0b0000; // Unknown
    }
  }

  // Encoding tables as constexpr for compile-time optimization
  static constexpr char decode_genome(uint8_t bit) {
    switch (bit) {
    case 0b0001:
      return 'A';
    case 0b0010:
      return 'C';
    case 0b0100:
      return 'G';
    case 0b1000:
      return 'T';
    case 0b0101:
      return 'R'; // A or G
    case 0b1010:
      return 'Y'; // C or T
    case 0b0110:
      return 'S'; // G or C
    case 0b1001:
      return 'W'; // A or T
    case 0b1100:
      return 'K'; // G or T
    case 0b0011:
      return 'M'; // A or C
    case 0b1110:
      return 'B'; // C or G or T
    case 0b1101:
      return 'D'; // A or G or T
    case 0b1011:
      return 'H'; // A or C or T
    case 0b0111:
      return 'V'; // A or C or G
    default:
      return '0'; // Unknown
    }
  }

  static constexpr uint8_t encode_pam(char c) {
    switch (c) {
    case 'A':
      return 0b0001;
    case 'C':
      return 0b0010;
    case 'G':
      return 0b0100;
    case 'T':
      return 0b1000;
    case 'N':
      return 0b1111; // Any (wildcard)
    case 'R':
      return 0b0101;
    case 'Y':
      return 0b1010;
    case 'S':
      return 0b0110;
    case 'W':
      return 0b1001;
    case 'K':
      return 0b1100;
    case 'M':
      return 0b0011;
    case 'B':
      return 0b1110;
    case 'D':
      return 0b1101;
    case 'H':
      return 0b1011;
    case 'V':
      return 0b0111;
    default:
      return 0b0000;
    }
  }

  static constexpr char decode_pam(uint8_t bit) {
    switch (bit) {
    case 0b0001:
      return 'A';
    case 0b0010:
      return 'C';
    case 0b0100:
      return 'G';
    case 0b1000:
      return 'T';
    case 0b1111:
      return 'N'; // Any (wildcard)
    case 0b0101:
      return 'R'; // A or G
    case 0b1010:
      return 'Y'; // C or T
    case 0b0110:
      return 'S'; // G or C
    case 0b1001:
      return 'W'; // A or T
    case 0b1100:
      return 'K'; // G or T
    case 0b0011:
      return 'M'; // A or C
    case 0b1110:
      return 'B'; // C or G or T
    case 0b1101:
      return 'D'; // A or G or T
    case 0b1011:
      return 'H'; // A or C or T
    case 0b0111:
      return 'V'; // A or C or G
    default:
      return '\0';
    }
  }

  static constexpr char complement(char c) {
    switch (c) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'U':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'R':
      return 'Y';
    case 'Y':
      return 'R';
    case 'M':
      return 'K';
    case 'K':
      return 'M';
    case 'H':
      return 'D';
    case 'D':
      return 'H';
    case 'B':
      return 'V';
    case 'V':
      return 'B';
    case 'S':
      return 'S';
    case 'W':
      return 'W';
    default:
      return c;
    }
  }
};

// Reverse complement a sequence
std::string reverse_complement(std::string_view seq);

} // namespace pam

#endif // NUCLEOTIDE_ENCODING_HPP