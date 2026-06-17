#include "search_configuration.hpp"

#include <stdexcept>
#include <string>

namespace crispritz {

// =========================================================================
// OutputFormat free functions
// =========================================================================

std::string_view to_string(OutputFormat fmt) noexcept {
  switch (fmt) {
  case OutputFormat::Tsv:
    return "tsv";
  case OutputFormat::Targets:
    return "targets";
  }
  return "tsv"; // unreachable for valid enum values
}

OutputFormat output_format_from_string(std::string_view name) {
  if (name == "tsv")
    return OutputFormat::Tsv;
  if (name == "targets")
    return OutputFormat::Targets;
  throw std::invalid_argument(
      "output_format_from_string: expected \"tsv\" or \"targets\", got \"" +
      std::string(name) + '"');
}

// =========================================================================
// OutputMode free functions
// =========================================================================

std::string_view to_string(OutputMode mode) noexcept {
  switch (mode) {
  case OutputMode::TargetsOnly:
    return "targets";
  case OutputMode::ProfileOnly:
    return "profile";
  case OutputMode::Both:
    return "both";
  }
  return "both"; // unreachable for valid enum values
}

OutputMode output_mode_from_string(std::string_view name) {
  if (name == "targets")
    return OutputMode::TargetsOnly;
  if (name == "profile")
    return OutputMode::ProfileOnly;
  if (name == "both")
    return OutputMode::Both;
  throw std::invalid_argument(
      "output_mode_from_string: expected \"targets\", \"profile\", or "
      "\"both\", got \"" +
      std::string(name) + '"');
}

// =========================================================================
// SearchConfiguration
// =========================================================================

SearchConfiguration::SearchConfiguration(int max_mismatches, int max_bulges_dna,
                                         int max_bulges_rna, int threads,
                                         OutputFormat output_format,
                                         OutputMode output_mode)
    : max_mismatches_(max_mismatches), max_bulges_dna_(max_bulges_dna),
      max_bulges_rna_(max_bulges_rna), threads_(threads),
      output_format_(output_format), output_mode_(output_mode) {
  // Only standalone-checkable invariants are validated here. The
  // edit-budget-vs-guide-length relationship is deferred to the index
  // loader, where guide_length becomes known (see header documentation).

  if (max_mismatches_ < 0)
    throw std::invalid_argument(
        "SearchConfiguration: max_mismatches must be >= 0, got " +
        std::to_string(max_mismatches_));

  if (max_bulges_dna_ < 0)
    throw std::invalid_argument(
        "SearchConfiguration: max_bulges_dna must be >= 0, got " +
        std::to_string(max_bulges_dna_));

  if (max_bulges_rna_ < 0)
    throw std::invalid_argument(
        "SearchConfiguration: max_bulges_rna must be >= 0, got " +
        std::to_string(max_bulges_rna_));

  if (threads_ < 1)
    throw std::invalid_argument(
        "SearchConfiguration: threads must be >= 1, got " +
        std::to_string(threads_));
}

} // namespace crispritz