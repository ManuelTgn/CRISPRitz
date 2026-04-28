/**
 * @file tst_loader.cpp
 * @brief Implementation of TSTLoader – binary TST file deserialization.
 */

#include "tst_loader.hpp"

#include <cassert>
#include <cstring>
#include <stdexcept>

namespace crispritz
{

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TSTLoader::TSTLoader(std::string path, int pam_size)
    : path_(std::move(path)), pam_size_(pam_size)
{
}

// ---------------------------------------------------------------------------
// Static helpers
// ---------------------------------------------------------------------------

char TSTLoader::nibble_to_char(uint8_t nibble) noexcept
{
    switch (nibble & 0x0F)
    {
    case 0x1: return 'A';
    case 0x2: return 'C';
    case 0x4: return 'G';
    case 0x8: return 'T';
    case 0x5: return 'R';
    case 0xA: return 'Y';
    case 0x6: return 'S';
    case 0x9: return 'W';
    case 0xC: return 'K';
    case 0x3: return 'M';
    case 0xE: return 'B';
    case 0xD: return 'D';
    case 0xB: return 'H';
    case 0x7: return 'V';
    case 0xF: return '_'; // sentinel
    default:  return '0'; // null child
    }
}

char TSTLoader::pam_nibble_to_char(uint8_t nibble) noexcept
{
    // High nibble of each byte encodes the PAM nucleotide
    switch (nibble & 0x0F)
    {
    case 0x1: return 'A';
    case 0x2: return 'C';
    case 0x4: return 'G';
    case 0x8: return 'T';
    case 0x5: return 'R';
    case 0xA: return 'Y';
    case 0x6: return 'S';
    case 0x9: return 'W';
    case 0xC: return 'K';
    case 0x3: return 'M';
    case 0xE: return 'B';
    case 0xD: return 'D';
    case 0xB: return 'H';
    case 0x7: return 'V';
    default:  return '\0';
    }
}

void TSTLoader::unpack_byte(char chars[2], std::bitset<4> bits[2], char raw) noexcept
{
    const auto byte = static_cast<unsigned char>(raw);

    // Sentinel: 0xF0 in the high nibble means '_' (leaf pointer marker)
    if ((byte & 0xF0) == 0xF0)
    {
        chars[0] = '_';
        bits[0]  = std::bitset<4>(0xF);
        // Low nibble is meaningless when high nibble is sentinel
        chars[1] = '0';
        bits[1]  = std::bitset<4>(0);
        return;
    }

    // High nibble → first character
    const uint8_t hi = (byte >> 4) & 0x0F;
    chars[0] = nibble_to_char(hi);
    bits[0]  = std::bitset<4>(hi);

    // Low nibble → second character
    const uint8_t lo = byte & 0x0F;
    chars[1] = nibble_to_char(lo);
    bits[1]  = std::bitset<4>(lo);
}

// ---------------------------------------------------------------------------
// Recursive deserialization
// ---------------------------------------------------------------------------

void TSTLoader::deserialize_node(std::ifstream& file, char chars[2],
                                  std::uint8_t bits[2], char& in,
                                  uint8_t& flag, int& node_idx)
{
    const int cur = node_idx;

    // Assign the current nibble to this node
    nodes_[cur].splitchar     = chars[flag];
    nodes_[cur].splitchar_enc = bits[flag];

    // Advance flag; read next byte if we consumed both nibbles
    auto advance = [&]() {
        if (flag == 1)
        {
            file.get(in);
            unpack_byte(chars, bits, in);
            flag = 0;
        }
        else
        {
            ++flag;
        }
    };

    advance();

    // lokid
    if (chars[flag] != '0')
    {
        nodes_[cur].lokid = ++node_idx;
        deserialize_node(file, chars, bits, in, flag, node_idx);
    }

    advance();

    // hikid
    if (chars[flag] != '0')
    {
        nodes_[cur].hikid = ++node_idx;
        deserialize_node(file, chars, bits, in, flag, node_idx);
    }

    advance();

    // eqkid: either a subtree or a leaf pointer (sentinel '_')
    if (chars[flag] == '_')
    {
        ++flag; // consume sentinel nibble
        int leaf_ptr = 0;
        file.read(reinterpret_cast<char*>(&leaf_ptr), sizeof(int));
        nodes_[cur].eqkid = leaf_ptr;
    }
    else
    {
        nodes_[cur].eqkid = ++node_idx;
        deserialize_node(file, chars, bits, in, flag, node_idx);
    }
}

// ---------------------------------------------------------------------------
// Public load()
// ---------------------------------------------------------------------------

void TSTLoader::load()
{
    std::ifstream file(path_, std::ios::in | std::ios::binary);
    if (!file.is_open())
        throw std::runtime_error("TSTLoader: cannot open index file: " + path_);

    // ---- header ------------------------------------------------------------
    int num_leaves = 0;
    file.read(reinterpret_cast<char*>(&num_leaves), sizeof(int));

    int stored_guide_len = 0;
    file.read(reinterpret_cast<char*>(&stored_guide_len), sizeof(int));
    // offset_guide_len_ is set by the caller after knowing the actual guide len;
    // here we store the raw value; the Searcher will compute the delta.
    offset_guide_len_ = stored_guide_len;

    // ---- leaf array --------------------------------------------------------
    leaves_.resize(static_cast<size_t>(num_leaves));

    for (int i = 0; i < num_leaves; ++i)
    {
        TSTLeaf& leaf = leaves_[i];

        file.read(reinterpret_cast<char*>(&leaf.guide_index), sizeof(int));

        // Decode packed PAM bytes: each byte holds two 4-bit IUPAC codes.
        // pam_size_ characters → ceil(pam_size_/2) bytes.
        leaf.pam_dna.resize(pam_size_);
        leaf.pam_dna_bit.resize(pam_size_);

        char in = 0;
        file.get(in);

        int k = 0; // tracks how many nibbles consumed from the current byte
        for (int j = pam_size_ - 1; j >= 0; --j)
        {
            if (k == 2)
            {
                file.get(in);
                k = 0;
            }

            // When PAM length is odd and we are at the last character,
            // shift to align the high nibble correctly.
            if (j == 0 && (pam_size_ % 2 != 0))
                in = static_cast<char>(static_cast<unsigned char>(in) << 4);

            const uint8_t mask = static_cast<unsigned char>(in) & 0xF0;
            // Shift next nibble into position
            in = static_cast<char>(static_cast<unsigned char>(in) << 4);

            char   pam_char = '\0';
            uint8_t pam_bits = 0;

            switch (mask)
            {
            case 0x10: pam_char = 'A'; pam_bits = 0b0001; break;
            case 0x20: pam_char = 'C'; pam_bits = 0b0010; break;
            case 0x40: pam_char = 'G'; pam_bits = 0b0100; break;
            case 0x80: pam_char = 'T'; pam_bits = 0b1000; break;
            case 0x50: pam_char = 'R'; pam_bits = 0b0101; break;
            case 0xA0: pam_char = 'Y'; pam_bits = 0b1010; break;
            case 0x60: pam_char = 'S'; pam_bits = 0b0110; break;
            case 0x90: pam_char = 'W'; pam_bits = 0b1001; break;
            case 0xC0: pam_char = 'K'; pam_bits = 0b1100; break;
            case 0x30: pam_char = 'M'; pam_bits = 0b0011; break;
            case 0xE0: pam_char = 'B'; pam_bits = 0b1110; break;
            case 0xD0: pam_char = 'D'; pam_bits = 0b1101; break;
            case 0xB0: pam_char = 'H'; pam_bits = 0b1011; break;
            case 0x70: pam_char = 'V'; pam_bits = 0b0111; break;
            default: break; // 0x00 → unknown, leave as '\0'
            }

            leaf.pam_dna[j]     = pam_char;
            leaf.pam_dna_bit[j] = std::bitset<4>(pam_bits);
            ++k;
        }

        // next-leaf chain link
        char chain_marker = 0;
        file.get(chain_marker);
        if (chain_marker == '0')
        {
            leaf.next = 0;
        }
        else
        {
            file.read(reinterpret_cast<char*>(&leaf.next), sizeof(int));
        }
    }

    // ---- node array --------------------------------------------------------
    int num_nodes = 0;
    file.read(reinterpret_cast<char*>(&num_nodes), sizeof(int));

    nodes_.resize(static_cast<size_t>(num_nodes));

    // Bootstrap the deserialization: read the first byte and set up state.
    char in = 0;
    file.get(in);

    char          chars[2];
    std::bitset<4> bits[2];
    unpack_byte(chars, bits, in);

    uint8_t flag     = 0; // start at the high nibble (index 0)
    int     node_idx = 0;

    deserialize_node(file, chars, bits, in, flag, node_idx);

    file.close();
}

} // namespace crispritz