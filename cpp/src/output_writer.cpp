#include "output_writer.hpp"

#include "offtarget.hpp"
#include "profile_data.hpp"

#include <fstream>
#include <iomanip> // std::setprecision, std::fixed
#include <memory>
#include <numeric> // std::accumulate
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace crispritz {

// =========================================================================
// TsvFormatter
// =========================================================================

std::string TsvFormatter::header() const {
  const std::vector<std::string> cols = OffTarget::tsv_header();
  std::string line;
  for (std::size_t i = 0; i < cols.size(); ++i) {
    if (i != 0)
      line += '\t';
    line += cols[i];
  }
  return line;
}

std::string TsvFormatter::format_row(const OffTarget &ot) const {
  return ot.to_tsv_row('\t');
}

// =========================================================================
// TargetsFormatter
// =========================================================================
//
// Legacy "targets" column order:
//   bulge_type, grna, target, chrom, pos, strand,
//   mismatches, bulge_size, total
//
// bulge_size = bulge_dna + bulge_rna
// total      = mismatches + bulge_size  (== total_edit_distance)
// =========================================================================

std::string TargetsFormatter::header() const {
  return "bulge_type\tgrna\ttarget\tchrom\tpos\tstrand\t"
         "mismatches\tbulge_size\ttotal";
}

std::string TargetsFormatter::format_row(const OffTarget &ot) const {
  const int bulge_size = ot.bulge_dna() + ot.bulge_rna();
  const int total = ot.total_edit_distance();

  std::string row;
  row.reserve(ot.grna().size() + ot.target().size() + 64u);

  row += ot.bulge_type();
  row += '\t';
  row += ot.grna();
  row += '\t';
  row += ot.target();
  row += '\t';
  row += ot.chrom();
  row += '\t';
  row += std::to_string(ot.pos());
  row += '\t';
  row += to_char(ot.strand());
  row += '\t';
  row += std::to_string(ot.mismatches());
  row += '\t';
  row += std::to_string(bulge_size);
  row += '\t';
  row += std::to_string(total);

  return row;
}

// =========================================================================
// Formatter factory
// =========================================================================

std::unique_ptr<OffTargetFormatter> make_formatter(OutputFormat fmt) {
  switch (fmt) {
  case OutputFormat::Tsv:
    return std::make_unique<TsvFormatter>();
  case OutputFormat::Targets:
    return std::make_unique<TargetsFormatter>();
  }
  throw std::invalid_argument(
      "make_formatter: unrecognised OutputFormat value " +
      std::to_string(static_cast<int>(fmt)));
}

// =========================================================================
// OutputWriter
// =========================================================================

OutputWriter::OutputWriter(OutputFormat fmt)
    : formatter_(make_formatter(fmt)) {}

OutputWriter::OutputWriter(std::unique_ptr<OffTargetFormatter> formatter)
    : formatter_(std::move(formatter)) {
  if (!formatter_)
    throw std::invalid_argument("OutputWriter: formatter must not be null");
}

namespace {
/**
 * @brief Throw if the stream is in a fail/bad state.
 * @param os    Stream to check.
 * @param where Context string for the error message.
 */
void check_stream(const std::ostream &os, const char *where) {
  if (!os.good())
    throw std::runtime_error(std::string("OutputWriter: write failed (") +
                             where + ')');
}
} // namespace

std::size_t OutputWriter::write(const std::vector<OffTarget> &records,
                                std::ostream &os, bool write_header) const {
  if (!os.good())
    throw std::runtime_error("OutputWriter: output stream is not writable");

  if (write_header) {
    os << formatter_->header() << '\n';
    check_stream(os, "header");
  }

  std::size_t written = 0;
  for (const OffTarget &ot : records) {
    os << formatter_->format_row(ot) << '\n';
    check_stream(os, "row");
    ++written;
  }
  return written;
}

std::size_t OutputWriter::write(const SearchResult &result, std::ostream &os,
                                bool write_header) const {
  if (!os.good())
    throw std::runtime_error("OutputWriter: output stream is not writable");

  if (write_header) {
    os << formatter_->header() << '\n';
    check_stream(os, "header");
  }

  std::size_t written = 0;
  for (const std::vector<OffTarget> &guide_hits : result.hits_by_guide) {
    for (const OffTarget &ot : guide_hits) {
      os << formatter_->format_row(ot) << '\n';
      check_stream(os, "row");
      ++written;
    }
  }
  return written;
}

std::size_t OutputWriter::write_to_file(const std::vector<OffTarget> &records,
                                        const std::string &path) const {
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out.is_open())
    throw std::runtime_error("OutputWriter: cannot open output file: " + path);

  return write(records, out, /*write_header=*/true);
}

std::size_t OutputWriter::write_to_file(const SearchResult &result,
                                        const std::string &path) const {
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out.is_open())
    throw std::runtime_error("OutputWriter: cannot open output file: " + path);

  return write(result, out, /*write_header=*/true);
}

// =========================================================================
// OutputWriter::Session — streaming, threshold-flushed writing
// =========================================================================

OutputWriter::Session::Session(const OffTargetFormatter &formatter,
                               std::ostream &os, std::size_t threshold,
                               bool write_header)
    : formatter_(formatter), owned_os_(nullptr), os_(os),
      threshold_(threshold) {
  if (threshold_ == 0)
    throw std::invalid_argument("OutputWriter::Session: threshold must be > 0");
  if (!os_.good())
    throw std::runtime_error(
        "OutputWriter::Session: output stream is not writable");

  buffer_.reserve(threshold_);

  if (write_header) {
    os_ << formatter_.header() << '\n';
    check_stream(os_, "header");
  }
}

OutputWriter::Session::Session(const OffTargetFormatter &formatter,
                               std::unique_ptr<std::ostream> owned_os,
                               std::size_t threshold, bool write_header)
    : formatter_(formatter), owned_os_(std::move(owned_os)), os_(*owned_os_),
      threshold_(threshold) {
  if (threshold_ == 0)
    throw std::invalid_argument("OutputWriter::Session: threshold must be > 0");
  if (!os_.good())
    throw std::runtime_error(
        "OutputWriter::Session: output stream is not writable");

  buffer_.reserve(threshold_);

  if (write_header) {
    os_ << formatter_.header() << '\n';
    check_stream(os_, "header");
  }
}

OutputWriter::Session::~Session() {
  // RAII safety net: flush whatever remains. A destructor must not throw,
  // so write errors here are suppressed. Call close() explicitly to have
  // such errors surface as exceptions.
  if (!closed_) {
    try {
      flush();
    } catch (...) { /* swallow — cannot propagate from a destructor */
    }
  }
}

void OutputWriter::Session::add(const OffTarget &ot) {
  buffer_.push_back(ot);
  if (buffer_.size() >= threshold_)
    flush();
}

void OutputWriter::Session::add_all(const std::vector<OffTarget> &records) {
  for (const OffTarget &ot : records)
    add(ot);
}

std::size_t OutputWriter::Session::flush() {
  if (buffer_.empty())
    return 0;

  for (const OffTarget &ot : buffer_) {
    os_ << formatter_.format_row(ot) << '\n';
    check_stream(os_, "row");
  }
  os_.flush();
  check_stream(os_, "flush");

  const std::size_t n = buffer_.size();
  written_ += n;
  buffer_.clear();
  return n;
}

std::size_t OutputWriter::Session::close() {
  if (!closed_) {
    flush();
    closed_ = true;
  }
  return written_;
}

OutputWriter::Session OutputWriter::open_session(std::ostream &os,
                                                 std::size_t threshold,
                                                 bool write_header) const {
  return Session(*formatter_, os, threshold, write_header);
}

OutputWriter::Session
OutputWriter::open_session_to_file(const std::string &path,
                                   std::size_t threshold) const {
  auto out =
      std::make_unique<std::ofstream>(path, std::ios::out | std::ios::trunc);
  if (!out->is_open())
    throw std::runtime_error("OutputWriter: cannot open output file: " + path);

  return Session(*formatter_, std::move(out), threshold, /*write_header=*/true);
}

// =========================================================================
// ProfileWriter — private helpers
// =========================================================================

std::string ProfileWriter::annotated_guide(const GuideProfile &p) {
  const std::string pam_ns(static_cast<std::size_t>(p.pam_len), 'N');
  return p.pam_at_start ? (pam_ns + p.guide) : (p.guide + pam_ns);
}

std::ofstream ProfileWriter::open_or_throw(const std::string &path) {
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out.is_open())
    throw std::runtime_error("ProfileWriter: cannot open output file: " + path);
  return out;
}

void ProfileWriter::check_stream(const std::ostream &os, const char *where) {
  if (!os.good())
    throw std::runtime_error(std::string("ProfileWriter: write failed (") +
                             where + ')');
}

void ProfileWriter::write_mm_header(int guide_len, int max_mm,
                                    std::ostream &os) {
  os << "GUIDE\t";
  for (int i = 0; i < guide_len; ++i)
    os << "BP\t";
  os << "\tONT\tOFFT\tMM/OFFT\t\t";
  for (int i = 0; i <= max_mm; ++i)
    os << i << "MM\t";
  os << '\n';
}

void ProfileWriter::write_bulge_header(int guide_len, int max_mm, int max_bulge,
                                       std::ostream &os) {
  os << "GUIDE\t";
  for (int i = 0; i < guide_len; ++i)
    os << "BP\t";
  os << "\tONT\tOFFT\tMM/OFFT\t\t";
  for (int mm = 0; mm <= max_mm; ++mm)
    for (int b = 1; b <= max_bulge; ++b)
      os << mm << "MM(" << b << ")\t";
  os << '\n';
}

// =========================================================================
// ProfileWriter — write_profile  (.profile.xls)
// =========================================================================

void ProfileWriter::write_profile(const std::vector<GuideProfile> &profiles,
                                  std::ostream &os) const {
  if (profiles.empty())
    return;

  const int guide_len = profiles.front().guide_len;
  const int max_mm = profiles.front().max_mm;

  write_mm_header(guide_len, max_mm, os);
  check_stream(os, "profile header");

  for (const GuideProfile &p : profiles) {
    // OFFT = all hits with ≥ 1 mismatch (MM-only channel).
    int offt_total = 0;
    int sum_mm = 0;
    for (int mm = 1; mm <= p.max_mm; ++mm) {
      offt_total += p.offt_by_mm[mm];
      sum_mm += p.offt_by_mm[mm] * mm;
    }
    const double mm_per_offt =
        (offt_total > 0) ? static_cast<double>(sum_mm) / offt_total : 0.0;

    os << annotated_guide(p) << '\t';
    for (int i = 0; i < p.guide_len; ++i)
      os << p.pos_mm_count[i] << '\t';
    os << '\t' << p.ont_count << '\t' << offt_total << '\t' << std::fixed
       << std::setprecision(6) << mm_per_offt << "\t\t";
    for (int mm = 0; mm <= p.max_mm; ++mm)
      os << p.offt_by_mm[mm] << '\t';
    os << '\n';
    check_stream(os, "profile row");
  }
}

// =========================================================================
// ProfileWriter — write_extended_profile  (.extended_profile.xls)
// =========================================================================

void ProfileWriter::write_extended_profile(
    const std::vector<GuideProfile> &profiles, std::ostream &os) const {
  static const char *NUC_LABELS[4] = {"NUCLEOTIDE A", "NUCLEOTIDE C",
                                      "NUCLEOTIDE G", "NUCLEOTIDE T"};

  for (const GuideProfile &p : profiles) {
    // FASTA-style guide header
    os << '>' << annotated_guide(p) << '\n';

    // Column-name row: leading tab + BP×L + TARGETS
    os << '\t';
    for (int i = 0; i < p.guide_len; ++i)
      os << "BP\t";
    os << "TARGETS\n";
    check_stream(os, "extended header");

    for (int mm = 0; mm <= p.max_mm; ++mm) {
      // Total-mismatches row
      os << mm << " MISMATCHES";
      for (int pos = 0; pos < p.guide_len; ++pos)
        os << '\t' << p.ext_total_by_mm[mm][pos];
      os << '\t' << p.offt_complete_by_mm[mm] << '\n';

      // Per-nucleotide rows (A, C, G, T)
      for (int nuc = 0; nuc < 4; ++nuc) {
        os << NUC_LABELS[nuc];
        for (int pos = 0; pos < p.guide_len; ++pos)
          os << '\t' << p.ext_mm_nuc_pos[mm][nuc][pos];
        os << '\n';
      }

      // DNA-bulge row
      os << "Bulge DNA";
      for (int pos = 0; pos < p.guide_len; ++pos)
        os << '\t' << p.ext_dna_by_mm_pos[mm][pos];
      os << '\n';

      // RNA-bulge row
      os << "Bulge RNA";
      for (int pos = 0; pos < p.guide_len; ++pos)
        os << '\t' << p.ext_rna_by_mm_pos[mm][pos];
      os << '\n';

      os << '\n'; // blank separator between mm blocks
      check_stream(os, "extended mm block");
    }
  }
}

// =========================================================================
// ProfileWriter — write_profile_dna  (.profile_dna.xls)
// =========================================================================

void ProfileWriter::write_profile_dna(const std::vector<GuideProfile> &profiles,
                                      std::ostream &os) const {
  if (profiles.empty())
    return;

  const int guide_len = profiles.front().guide_len;
  const int max_mm = profiles.front().max_mm;
  const int max_bulge = profiles.front().max_bulge_dna;

  // Skip header / data entirely when bulges were not searched.
  if (max_bulge <= 0)
    return;

  write_bulge_header(guide_len, max_mm, max_bulge, os);
  check_stream(os, "profile_dna header");

  for (const GuideProfile &p : profiles) {
    int offt_total = 0;
    int sum_mm = 0;
    for (int mm = 0; mm <= p.max_mm; ++mm) {
      for (int b = 0; b < p.max_bulge_dna; ++b) {
        offt_total += p.offt_dna[mm][b];
        sum_mm += p.offt_dna[mm][b] * mm;
      }
    }
    // Subtract 0MM (on-target) from offt_total.
    int ont = p.ont_count_dna;
    int offt = offt_total - ont;
    const double mm_per_offt =
        (offt > 0) ? static_cast<double>(sum_mm) / offt : 0.0;

    os << annotated_guide(p) << '\t';
    // BP cells: mm_in_dna(bulge_count) format
    for (int i = 0; i < p.guide_len; ++i)
      os << p.pos_mm_in_dna[i] << '(' << p.pos_bulge_dna[i] << ")\t";
    os << '\t' << ont << '\t' << offt << '\t' << std::fixed
       << std::setprecision(6) << mm_per_offt << "\t\t";
    for (int mm = 0; mm <= p.max_mm; ++mm)
      for (int b = 0; b < p.max_bulge_dna; ++b)
        os << p.offt_dna[mm][b] << '\t';
    os << '\n';
    check_stream(os, "profile_dna row");
  }
}

// =========================================================================
// ProfileWriter — write_profile_rna  (.profile_rna.xls)
// =========================================================================

void ProfileWriter::write_profile_rna(const std::vector<GuideProfile> &profiles,
                                      std::ostream &os) const {
  if (profiles.empty())
    return;

  const int guide_len = profiles.front().guide_len;
  const int max_mm = profiles.front().max_mm;
  const int max_bulge = profiles.front().max_bulge_rna;

  if (max_bulge <= 0)
    return;

  write_bulge_header(guide_len, max_mm, max_bulge, os);
  check_stream(os, "profile_rna header");

  for (const GuideProfile &p : profiles) {
    int offt_total = 0;
    int sum_mm = 0;
    for (int mm = 0; mm <= p.max_mm; ++mm) {
      for (int b = 0; b < p.max_bulge_rna; ++b) {
        offt_total += p.offt_rna[mm][b];
        sum_mm += p.offt_rna[mm][b] * mm;
      }
    }
    int ont = p.ont_count_rna;
    int offt = offt_total - ont;
    const double mm_per_offt =
        (offt > 0) ? static_cast<double>(sum_mm) / offt : 0.0;

    os << annotated_guide(p) << '\t';
    for (int i = 0; i < p.guide_len; ++i)
      os << p.pos_mm_in_rna[i] << '(' << p.pos_bulge_rna[i] << ")\t";
    os << '\t' << ont << '\t' << offt << '\t' << std::fixed
       << std::setprecision(6) << mm_per_offt << "\t\t";
    for (int mm = 0; mm <= p.max_mm; ++mm)
      for (int b = 0; b < p.max_bulge_rna; ++b)
        os << p.offt_rna[mm][b] << '\t';
    os << '\n';
    check_stream(os, "profile_rna row");
  }
}

// =========================================================================
// ProfileWriter — write_profile_complete  (.profile_complete.xls)
// =========================================================================

void ProfileWriter::write_profile_complete(
    const std::vector<GuideProfile> &profiles, std::ostream &os) const {
  if (profiles.empty())
    return;

  const int guide_len = profiles.front().guide_len;
  const int max_mm = profiles.front().max_mm;

  write_mm_header(guide_len, max_mm, os);
  check_stream(os, "profile_complete header");

  for (const GuideProfile &p : profiles) {
    int offt_total = 0;
    int sum_mm = 0;
    for (int mm = 1; mm <= p.max_mm; ++mm) {
      offt_total += p.offt_complete_by_mm[mm];
      sum_mm += p.offt_complete_by_mm[mm] * mm;
    }
    const double mm_per_offt =
        (offt_total > 0) ? static_cast<double>(sum_mm) / offt_total : 0.0;

    os << annotated_guide(p) << '\t';
    for (int i = 0; i < p.guide_len; ++i)
      os << p.pos_mm_complete[i] << '\t';
    os << '\t' << p.ont_count_complete << '\t' << offt_total << '\t'
       << std::fixed << std::setprecision(6) << mm_per_offt << "\t\t";
    for (int mm = 0; mm <= p.max_mm; ++mm)
      os << p.offt_complete_by_mm[mm] << '\t';
    os << '\n';
    check_stream(os, "profile_complete row");
  }
}

// =========================================================================
// ProfileWriter — write_all_profiles
// =========================================================================

void ProfileWriter::write_all_profiles(
    const std::vector<GuideProfile> &profiles,
    const std::string &path_stem) const {
  {
    auto f = open_or_throw(path_stem + ".profile.xls");
    write_profile(profiles, f);
  }
  {
    auto f = open_or_throw(path_stem + ".extended_profile.xls");
    write_extended_profile(profiles, f);
  }
  {
    auto f = open_or_throw(path_stem + ".profile_dna.xls");
    write_profile_dna(profiles, f);
  }
  {
    auto f = open_or_throw(path_stem + ".profile_rna.xls");
    write_profile_rna(profiles, f);
  }
  {
    auto f = open_or_throw(path_stem + ".profile_complete.xls");
    write_profile_complete(profiles, f);
  }
}

} // namespace crispritz