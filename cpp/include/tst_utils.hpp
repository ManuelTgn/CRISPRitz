#ifndef TST_UTILS_HPP
#define TST_UTILS_HPP

#include "nucleotide_encoding.hpp"

#include <cstdint>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace crispritz
{

    namespace iupac 
    {
        uint8_t encode_genome(char c) {
            return pam::NucleotideEncoder::encode_genome(c);
        }

        char decode_genome(uint8_t bit) {
            return pam::NucleotideEncoder::decode_genome(bit);
        }

        char complement(char c) {
            return pam::NucleotideEncoder::complement(c);
        }

        std::string reverse_complement(std::string_view seq) {
            return pam::reverse_complement(seq);
        }
    }

    /**
     * @brief Sentinel high-nibble value used to signal end-of-node in the binary
     *        format (equivalent to the legacy '_' character in writePair).
     */
    constexpr uint8_t SENTINEL_HIGH_NIBBLE = 0xF0;

    /**
     * @brief Sentinel low-nibble value used to signal end-of-node in the binary
     *        format (equivalent to the legacy '_' character in writePair).
     */
    constexpr uint8_t SENTINEL_LOW_NIBBLE = 0x0F;

    /**
     * @brief Sentinel byte written when a TST node has no child in that slot.
     *        Equivalent to the legacy '0' character written by serialize().
     */
    constexpr uint8_t NULL_CHILD_NIBBLE = 0x0;
    
    class Serializer
    {
      public:
        // Encoding tables as constexpr for compile-time optimization
        static constexpr uint8_t serialize_left_nuc(char c)
        {
            switch (c)
            {
                case 'A': return 0x10;
                case 'C': return 0x20;
                case 'G': return 0x40;
                case 'T': return 0x80;
                case '_': return SENTINEL_HIGH_NIBBLE;
                case 'R': return 0x50;
                case 'Y': return 0xA0;
                case 'S': return 0x60;
                case 'W': return 0x90;
                case 'K': return 0xC0;
                case 'M': return 0x30;
                case 'B': return 0xE0;
                case 'D': return 0xD0;
                case 'H': return 0xB0;
                case 'V': return 0x70;
                case '0': return NULL_CHILD_NIBBLE;
                default:  return NULL_CHILD_NIBBLE;
            }
        }

        static constexpr uint8_t serialize_right_nuc(char c)
        {
            switch (c)
            {
                case 'A': return 0x1; 
                case 'C': return 0x2;
                case 'G': return 0x4;
                case 'T': return 0x8;
                case '_': return SENTINEL_LOW_NIBBLE;
                case 'R': return 0x05;
                case 'Y': return 0x0A;
                case 'S': return 0x06;
                case 'W': return 0x09;
                case 'K': return 0x0C;
                case 'M': return 0x03;
                case 'B': return 0x0E;
                case 'D': return 0x0D;
                case 'H': return 0x0B;
                case 'V': return 0x07;
                case '0': return NULL_CHILD_NIBBLE;
                default:  return NULL_CHILD_NIBBLE;
            }
        }

        static constexpr uint8_t serialize_pam(char c)
        {
            switch (c)
            {
                case 'A': return 0x1; 
                case 'C': return 0x2;
                case 'G': return 0x4;
                case 'T': return 0x8;
                case 'R': return 0x05;
                case 'Y': return 0x0A;
                case 'S': return 0x06;
                case 'W': return 0x09;
                case 'K': return 0x0C;
                case 'M': return 0x03;
                case 'B': return 0x0E;
                case 'D': return 0x0D;
                case 'H': return 0x0B;
                case 'V': return 0x07;
                default:  return NULL_CHILD_NIBBLE;
            }
        }
    };

    inline bool is_valid_genomic_window(std::string window) noexcept {
        for (char c : window) {
            // skip window containing 'N' or null elements
            if (c == 'N' || c == '\0' || c == '-')
                return false;
        }
        return true;
    }

    std::filesystem::path format_partition_filename(std::string outdir, std::string pam_seq, std::string chr_name, int part) {
        return std::filesystem::path(outdir) / (pam_seq + "_" + chr_name + "_" + std::to_string(part) + ".bin");
    }
    

} // namespace crispritz

#endif // TST_UTILS_HPP
