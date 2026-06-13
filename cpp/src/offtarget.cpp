#include "offtarget.hpp"

#include <stdexcept>
#include <string>

namespace crispritz {

// =========================================================================
// Free functions
// =========================================================================

Strand strand_from_char(char c) {
  if (c == '+')
    return Strand::Forward;
  if (c == '-')
    return Strand::Reverse;
  throw std::invalid_argument(
      std::string("strand_from_char: expected '+' or '-', got '") + c + '\'');
}

// =========================================================================
// Constructor
// =========================================================================

OffTarget::OffTarget(std::string chrom, int32_t pos, Strand strand,
                     std::string grna, std::string target, int mismatches,
                     int bulge_dna, int bulge_rna)
    : chrom_(std::move(chrom)), pos_(pos), strand_(strand),
      grna_(std::move(grna)), target_(std::move(target)),
      mismatches_(mismatches), bulge_dna_(bulge_dna), bulge_rna_(bulge_rna) {
  // Validate after moving so error messages can still reference the values.
  if (chrom_.empty())
    throw std::invalid_argument("OffTarget: chrom must not be empty");

  if (pos_ <= 0)
    throw std::invalid_argument(
        "OffTarget: pos must be > 0 (1-based coordinate), got " +
        std::to_string(pos));

  if (grna_.empty())
    throw std::invalid_argument("OffTarget: grna must not be empty");

  if (target_.empty())
    throw std::invalid_argument("OffTarget: target must not be empty");

  if (mismatches_ < 0)
    throw std::invalid_argument("OffTarget: mismatches must be >= 0, got " +
                                std::to_string(mismatches));

  if (bulge_dna_ < 0)
    throw std::invalid_argument("OffTarget: bulge_dna must be >= 0, got " +
                                std::to_string(bulge_dna));

  if (bulge_rna_ < 0)
    throw std::invalid_argument("OffTarget: bulge_rna must be >= 0, got " +
                                std::to_string(bulge_rna));
}

// =========================================================================
// Derived properties
// =========================================================================

std::string OffTarget::bulge_type() const {
  // Derived from bulge_dna_ / bulge_rna_ — no separate stored field
  // to keep the class invariant simple (the three fields always agree).
  if (bulge_dna_ == 0 && bulge_rna_ == 0)
    return "X";
  if (bulge_dna_ > 0 && bulge_rna_ == 0)
    return "DNA";
  if (bulge_dna_ == 0 && bulge_rna_ > 0)
    return "RNA";
  return "DNA,RNA";
}

// =========================================================================
// Comparison operators
// =========================================================================

bool OffTarget::operator==(const OffTarget &other) const noexcept {
  // Locus identity only — guide sequence and edit counts are excluded.
  // See class-level documentation for the rationale.
  return chrom_ == other.chrom_ && pos_ == other.pos_ &&
         strand_ == other.strand_ && target_ == other.target_;
}

bool OffTarget::operator!=(const OffTarget &other) const noexcept {
  return !(*this == other);
}

bool OffTarget::operator<(const OffTarget &other) const noexcept {
  // Priority: chrom → pos → strand → total edit distance.
  //
  // std::tie cannot be used here because total_edit_distance() is a
  // computed value (not a stored l-value), so a cascaded if-chain is
  // clearer and avoids temporary storage.
  if (chrom_ != other.chrom_)
    return chrom_ < other.chrom_;

  if (pos_ != other.pos_)
    return pos_ < other.pos_;

  // Strand: '+' (ASCII 43) < '-' (ASCII 45).
  if (strand_ != other.strand_)
    return to_char(strand_) < to_char(other.strand_);

  return total_edit_distance() < other.total_edit_distance();
}

// =========================================================================
// Serialization
// =========================================================================

std::string OffTarget::to_tsv_row(char sep) const {
  // Build via string concatenation with a pre-reserved capacity.
  // Avoids the construction overhead and locale state of std::ostringstream
  // for what is a simple, fixed-schema row.
  std::string row;
  row.reserve(chrom_.size() + grna_.size() + target_.size() + 64u);

  row += chrom_;
  row += sep;
  row += std::to_string(pos_);
  row += sep;
  row += to_char(strand_);
  row += sep;
  row += grna_;
  row += sep;
  row += target_;
  row += sep;
  row += std::to_string(mismatches_);
  row += sep;
  row += std::to_string(bulge_dna_);
  row += sep;
  row += std::to_string(bulge_rna_);
  row += sep;
  row += bulge_type();

  return row;
}

std::vector<std::string> OffTarget::tsv_header() {
  // Column order mirrors to_tsv_row().  CFD score is deliberately absent;
  // the Python scoring layer appends it after computing it over the rows
  // returned from the C++ search.
  return {"chrom",      "pos",       "strand",    "grna",      "target",
          "mismatches", "bulge_dna", "bulge_rna", "bulge_type"};
}

// =========================================================================
// Debug / diagnostic
// =========================================================================

std::string OffTarget::locus_string() const {
  // E.g. "chr1:123456(+)"
  return chrom_ + ':' + std::to_string(pos_) + '(' + to_char(strand_) + ')';
}

} // namespace crispritz