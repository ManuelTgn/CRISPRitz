#ifndef TST_UTILS_HPP
#define TST_UTILS_HPP

#include "nucleotide_encoding.hpp" // NucleotideEncoder, pam::reverse_complement

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace crispritz {

/**
 * @brief 4-bit IUPAC nucleotide encoding used throughout TST construction.
 *
 * Bit layout (LSB first): bit0=A, bit1=C, bit2=G, bit3=T.
 * Ambiguous codes set the union of their constituent bases.
 * The sentinel value 0x00 represents a gap character '-' (used in bulge
 * alignment). The value 0x0F represents 'N' (any base) in PAM contexts and
 * also the '_' end-of-sequence marker in the serialized binary format.
 *
 * These values must stay stable: the binary .bin format written by saveTST
 * and read by searchOnTST encodes two nucleotides per byte using these exact
 * 4-bit patterns.
 */
namespace iupac {
// aliases into the canonical implementation so tst.cpp call sites
// do not need to change namespace, and there is exactly one encoding table
using pam::NucleotideEncoder;

constexpr uint8_t encode_genome(char c) noexcept {
  return NucleotideEncoder::encode_genome(c);
}

constexpr uint8_t encode_pam(char c) noexcept {
  return NucleotideEncoder::encode_genome(c);
}

constexpr char complement(char c) noexcept {
  return NucleotideEncoder::complement(c);
}

constexpr bool matches(uint8_t a, uint8_t b) noexcept { return (a & b) != 0; }
} // namespace iupac

// alias the already-implemented reverse_complement from the pam namespace
using pam::reverse_complement;

// ---------------------------------------------------------------------------
// Bit-packing helpers (two IUPAC nibbles per byte)
//
// The binary .bin format stores two 4-bit IUPAC encoded nucleotides per byte:
//   high nibble (bits 7-4) = first nucleotide
//   low  nibble (bits 3-0) = second nucleotide
//
// The value 0b1111 in the high nibble is the special end-of-sequence sentinel
// '_' used during TST serialization (writePair equivalent).
//
// These helpers centralize the packing/unpacking so that both the TST writer
// (tst.cpp) and any future reader share identical logic.
// ---------------------------------------------------------------------------

/**
 * @brief Pack two 4-bit IUPAC encodings into one byte.
 *
 * @param high  Encoding for the high nibble (first nucleotide).
 * @param low   Encoding for the low nibble  (second nucleotide).
 * @return      Packed byte: (high << 4) | low.
 */
constexpr uint8_t pack_nibbles(uint8_t high, uint8_t low) noexcept {
  return static_cast<uint8_t>((high << 4) | (low & 0b1111));
}

/**
 * @brief Extract the high nibble (first nucleotide) from a packed byte.
 * @param byte  Packed byte produced by pack_nibbles.
 * @return      4-bit IUPAC encoding of the first nucleotide.
 */
constexpr uint8_t high_nibble(uint8_t byte) noexcept {
  return static_cast<uint8_t>((byte >> 4) & 0b1111);
}

/**
 * @brief Extract the low nibble (second nucleotide) from a packed byte.
 * @param byte  Packed byte produced by pack_nibbles.
 * @return      4-bit IUPAC encoding of the second nucleotide.
 */
constexpr uint8_t low_nibble(uint8_t byte) noexcept {
  return static_cast<uint8_t>(byte & 0b1111);
}

/**
 * @brief Sentinel high-nibble value used to signal end-of-node in the binary
 *        format (equivalent to the legacy '_' character in writePair).
 */
constexpr uint8_t SENTINEL_NIBBLE = 0b1111;

/**
 * @brief Sentinel byte written when a TST node has no child in that slot.
 *        Equivalent to the legacy '0' character written by serialize().
 */
constexpr uint8_t NULL_CHILD_NIBBLE = 0b0000;

} // namespace crispritz

#endif // TST_UTILS_HPP
