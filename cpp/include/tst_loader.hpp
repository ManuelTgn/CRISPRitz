#pragma once

/**
 * @file tst_loader.hpp
 * @brief Deserializes a binary TST index file produced by the genome-indexing
 *        step into an in-memory node array and leaf array.
 *
 * @code
 * [int32 num_leaves]
 * [int32 guide_length]
 * for each leaf:
 *   [int32  guide_index]
 *   [ceil(pam_size/2) bytes  packed IUPAC PAM]
 *   ['0' | ('_' [int32 next])]
 * [int32 num_nodes]
 * [packed node stream -- pairs of 4-bit IUPAC chars]
 * @endcode
 */

#include "offtarget.hpp"
#include "tst.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace crispritz
{

/**
 * @brief Loads a single TST partition (.bin file) into memory.
 *
 * After construction call @c load(); the node array and leaf array can then
 * be accessed via @c nodes() and @c leaves().
 *
 * @par Thread Safety
 * Each thread should construct its own @c TSTLoader instance.
 */
class TSTLoader
{
  public:
    /**
     * @brief Constructs a loader for the given binary file path.
     * @param path   Absolute or relative path to the @c .bin index file.
     * @param pam_size  Length (in characters) of the PAM string stored in
     *                  each leaf (e.g. 3 for NGG).
     */
    explicit TSTLoader(std::string path, int pam_size);

    /**
     * @brief Parses the binary file and populates the node/leaf arrays.
     * @throws std::runtime_error if the file cannot be opened or is malformed.
     */
    void load();

    /**
     * @brief Returns the guide length delta stored in the index file.
     *
     * The legacy format stores the guide length that was used at build time.
     * Searching with a guide of different length requires an offset correction.
     *
     * @return guide_length_stored - current guide_length (may be zero)
     */
    int guide_length_offset() const noexcept { return offset_guide_len_; }

    /** @return Read-only reference to the deserialized node array. */
    const std::vector<TSTNode>& nodes() const noexcept { return nodes_; }

    /** @return Read-only reference to the deserialized leaf array. */
    const std::vector<TSTLeaf>& leaves() const noexcept { return leaves_; }

    /** @return Number of nodes in the deserialized tree. */
    int num_nodes() const noexcept { return static_cast<int>(nodes_.size()); }

    /** @return Number of leaves in the deserialized tree. */
    int num_leaves() const noexcept { return static_cast<int>(leaves_.size()); }

  private:
    // ---- helpers -----------------------------------------------------------

    /**
     * @brief Reads one byte from @c file_ and unpacks it into two IUPAC chars
     *        and their 4-bit encodings.
     *
     * The binary format packs pairs of nibbles: high 4 bits = first char,
     * low 4 bits = second char.  The sentinel value 0xF0 signals the '_'
     * (leaf-pointer) marker.
     *
     * @param[out] chars   Two-char buffer to receive the decoded characters.
     * @param[out] bits    Two-element bitset buffer to receive the encodings.
     * @param[in]  byte    The raw byte already read from the file.
     */
    static void unpack_byte(char chars[2], std::bitset<4> bits[2], char byte) noexcept;

    /**
     * @brief Recursively rebuilds the TST from the serialized nibble stream.
     *
     * Mirrors the legacy @c deSerialize function but operates on class members
     * instead of global state.
     *
     * @param[in,out] file      Open input file stream.
     * @param[in,out] chars     Two-char decode buffer shared across recursion.
     * @param[in,out] bits      Two-bitset decode buffer shared across recursion.
     * @param[in,out] in        Most recently read raw byte.
     * @param[in,out] flag      Current nibble index (0 or 1).
     * @param[in,out] node_idx  Index of the node being populated (auto-incremented).
     */
    void deserialize_node(std::ifstream& file, char chars[2], std::bitset<4> bits[2],
                          char& in, uint8_t& flag, int& node_idx);

    // ---- data members ------------------------------------------------------
    std::string           path_;
    int                   pam_size_;
    int                   offset_guide_len_{0};
    std::vector<TSTNode>  nodes_;
    std::vector<TSTLeaf>  leaves_;
};

} // namespace crispritz