#pragma once

#include "offtarget.hpp"            // crispritz::OffTarget
#include "profile_data.hpp"         // crispritz::GuideProfile
#include "search_configuration.hpp" // crispritz::OutputFormat
#include "tst_search.hpp"           // crispritz::SearchResult

#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace crispritz {

// =========================================================================
// OffTargetFormatter — abstract formatting strategy
// =========================================================================

/**
 * @brief Strategy interface that turns OffTarget records into text lines.
 *
 * A formatter encapsulates exactly one output layout: the header line (if
 * any) and the per-record row. It performs no I/O — it only produces
 * strings. This keeps the *what to write* (formatter) cleanly separated
 * from the *where/how to write it* (OutputWriter), and from the *what was
 * found* (the search layer).
 *
 * Implementations must be stateless and reentrant so a single formatter can
 * format an arbitrary number of records.
 */
class OffTargetFormatter {
public:
  virtual ~OffTargetFormatter() = default;

  /**
   * @brief The header line for this format, without a trailing newline.
   * @return Header string; empty if the format has no header.
   * @complexity O(number of columns).
   */
  [[nodiscard]] virtual std::string header() const = 0;

  /**
   * @brief Format one off-target as a single line, without a newline.
   * @param ot  The record to format.
   * @return    One formatted row.
   * @complexity O(length of the record's sequences).
   */
  [[nodiscard]] virtual std::string format_row(const OffTarget &ot) const = 0;

  /**
   * @brief Canonical name of the format ("tsv", "targets", ...).
   * @return A view into static storage.
   */
  [[nodiscard]] virtual std::string_view name() const noexcept = 0;
};

// =========================================================================
// Concrete formatters
// =========================================================================

/**
 * @brief Tab-separated values matching OffTarget::tsv_header() /
 *        OffTarget::to_tsv_row().
 *
 * Column order (9 columns):
 *   chrom, pos, strand, grna, target, mismatches,
 *   bulge_dna, bulge_rna, bulge_type
 *
 * This is the canonical search-side schema. The Python layer later maps
 * @c target → @c spacer and appends the CFD score column; the C++ writer is
 * deliberately not responsible for the scored final TSV.
 */
class TsvFormatter final : public OffTargetFormatter {
public:
  [[nodiscard]] std::string header() const override;
  [[nodiscard]] std::string format_row(const OffTarget &ot) const override;
  [[nodiscard]] std::string_view name() const noexcept override {
    return "tsv";
  }
};

/**
 * @brief Legacy CRISPRitz "targets" layout.
 *
 * The legacy detailedOutput "targets" file uses the same underlying fields
 * as the TSV but in the historical column order and with the bulge-type
 * column first, matching what downstream legacy tooling expects:
 *
 *   bulge_type, grna, target, chrom, pos, strand,
 *   mismatches, bulge_size, total
 *
 * where @c bulge_size = bulge_dna + bulge_rna and
 *       @c total       = mismatches + bulge_size.
 *
 * Field values themselves are identical to the TSV; only the column
 * selection and order differ.
 */
class TargetsFormatter final : public OffTargetFormatter {
public:
  [[nodiscard]] std::string header() const override;
  [[nodiscard]] std::string format_row(const OffTarget &ot) const override;
  [[nodiscard]] std::string_view name() const noexcept override {
    return "targets";
  }
};

/**
 * @brief Construct the formatter for a given OutputFormat.
 *
 * Single point of format dispatch. Adding a new format requires a new
 * OffTargetFormatter subclass and one new case here — nothing in
 * OutputWriter changes.
 *
 * @param fmt  The desired output format.
 * @return     An owning pointer to a formatter for @p fmt.
 * @throws std::invalid_argument if @p fmt is unrecognised.
 */
[[nodiscard]] std::unique_ptr<OffTargetFormatter>
make_formatter(OutputFormat fmt);

// =========================================================================
// OutputWriter
// =========================================================================

/**
 * @brief Writes off-target records to a stream or file in a chosen format.
 *
 * OutputWriter owns a formatter and handles all mechanical output concerns:
 * writing the header once, emitting one line per record, and (for the file
 * overload) opening and RAII-closing the destination. It is the single
 * component permitted to perform output formatting and file I/O; the search
 * layer knows nothing about either.
 *
 * Typical use:
 * @code
 *   OutputWriter writer(OutputFormat::Tsv);
 *   writer.write_to_file(result, "chr1.targets.tsv");
 * @endcode
 */
class OutputWriter {
public:
  /**
   * @brief Construct a writer for the given format.
   * @param fmt  Output format (default: OutputFormat::Tsv).
   * @throws std::invalid_argument if @p fmt is unrecognised.
   */
  explicit OutputWriter(OutputFormat fmt = OutputFormat::Tsv);

  /**
   * @brief Construct a writer with an explicit formatter (for testing or
   *        custom formats not in the OutputFormat enum).
   * @param formatter  Owning formatter pointer (must not be null).
   * @throws std::invalid_argument if @p formatter is null.
   */
  explicit OutputWriter(std::unique_ptr<OffTargetFormatter> formatter);

  // ---- stream API (testable without the filesystem) ------------------

  /**
   * @brief Write a flat list of records (header + one row each) to a stream.
   *
   * @param records  Records to write, in the order given.
   * @param os       Destination stream (must be good() on entry).
   * @param write_header  Emit the header line first (default true).
   * @return         The number of record rows written.
   * @throws std::runtime_error if the stream enters a fail state mid-write.
   */
  std::size_t write(const std::vector<OffTarget> &records, std::ostream &os,
                    bool write_header = true) const;

  /**
   * @brief Write a SearchResult (all guides' hits) to a stream.
   *
   * Records are written guide by guide, in guide order. A single header
   * is emitted once at the top.
   *
   * @param result       The search result to serialize.
   * @param os           Destination stream.
   * @param write_header Emit the header line first (default true).
   * @return             Total number of record rows written.
   * @throws std::runtime_error if the stream enters a fail state mid-write.
   */
  std::size_t write(const SearchResult &result, std::ostream &os,
                    bool write_header = true) const;

  // ---- file API ------------------------------------------------------

  /**
   * @brief Write a flat list of records to a file path.
   *
   * Opens @p path for writing (truncating), writes header + rows, and
   * closes via RAII. The directory must already exist.
   *
   * @param records  Records to write.
   * @param path     Destination file path.
   * @return         The number of record rows written.
   * @throws std::runtime_error if the file cannot be opened or a write fails.
   */
  std::size_t write_to_file(const std::vector<OffTarget> &records,
                            const std::string &path) const;

  /**
   * @brief Write a SearchResult to a file path.
   * @param result  The search result to serialize.
   * @param path    Destination file path.
   * @return        Total number of record rows written.
   * @throws std::runtime_error if the file cannot be opened or a write fails.
   */
  std::size_t write_to_file(const SearchResult &result,
                            const std::string &path) const;

  /** @return The canonical name of the active format. */
  [[nodiscard]] std::string_view format_name() const noexcept {
    return formatter_->name();
  }

  // ---- streaming / threshold-flush API -------------------------------

  /**
   * @brief Default number of buffered records before an automatic flush.
   *
   * Bounds peak memory during genome-wide searches: instead of holding an
   * entire partition's hits (potentially millions) in memory before
   * writing, a streaming Session accumulates at most this many records and
   * flushes them to the output stream when the buffer fills. This mirrors
   * the legacy "write after N targets" behaviour.
   */
  static constexpr std::size_t DEFAULT_FLUSH_THRESHOLD = 1'000'000;

  /**
   * @brief A streaming write session that flushes records in capped batches.
   *
   * Created via OutputWriter::open_session() / open_session_to_file(). The
   * session writes the header once on construction, buffers records via
   * add()/add_all(), and flushes automatically whenever the buffer reaches
   * the configured threshold. Any remaining buffered records are flushed
   * by close() or by the destructor (RAII), so no records are lost even on
   * early scope exit.
   *
   * The session holds at most @c threshold records in memory at once,
   * independent of how many records are ultimately written — this is what
   * makes genome-wide output memory-safe.
   *
   * A Session is move-only (it owns a buffer and references a stream); it
   * must not outlive the stream it writes to.
   */
  class Session {
  public:
    /**
     * @brief Begin a streaming session on @p os.
     * @param formatter    Formatter to use (borrowed; must outlive Session).
     * @param os           Destination stream (borrowed; must outlive Session).
     * @param threshold    Records buffered before an auto-flush (> 0).
     * @param write_header Emit the header immediately (default true).
     * @throws std::runtime_error     if the stream is not writable.
     * @throws std::invalid_argument  if threshold == 0.
     */
    Session(const OffTargetFormatter &formatter, std::ostream &os,
            std::size_t threshold, bool write_header = true);

    /**
     * @brief Optional owning variant: the session owns the file stream.
     *
     * Used by open_session_to_file() so the file is closed (after a final
     * flush) when the session is destroyed.
     */
    Session(const OffTargetFormatter &formatter,
            std::unique_ptr<std::ostream> owned_os, std::size_t threshold,
            bool write_header = true);

    ~Session();

    Session(Session &&) noexcept = default;
    Session &operator=(Session &&) noexcept = default;
    Session(const Session &) = delete;
    Session &operator=(const Session &) = delete;

    /**
     * @brief Buffer one record, auto-flushing if the threshold is hit.
     * @param ot  Record to add.
     * @throws std::runtime_error if an auto-flush write fails.
     */
    void add(const OffTarget &ot);

    /**
     * @brief Buffer many records (each may trigger an auto-flush).
     * @param records  Records to add, in order.
     */
    void add_all(const std::vector<OffTarget> &records);

    /**
     * @brief Flush all currently buffered records to the stream now.
     * @return Number of records flushed by this call.
     * @throws std::runtime_error if the stream enters a fail state.
     */
    std::size_t flush();

    /**
     * @brief Flush remaining records and mark the session complete.
     *
     * Idempotent. Called automatically by the destructor; call it
     * explicitly when you want write errors to surface as exceptions
     * (a destructor must not throw, so errors during destruction are
     * suppressed).
     *
     * @return Total records written over the session's lifetime.
     * @throws std::runtime_error if the final flush fails.
     */
    std::size_t close();

    /** @return Total records written so far (flushed). */
    [[nodiscard]] std::size_t total_written() const noexcept {
      return written_;
    }

    /** @return Records currently buffered, not yet flushed. */
    [[nodiscard]] std::size_t buffered() const noexcept {
      return buffer_.size();
    }

  private:
    const OffTargetFormatter &formatter_;
    std::unique_ptr<std::ostream>
        owned_os_; // non-null iff session owns the stream
    std::ostream &os_;
    std::size_t threshold_;
    std::vector<OffTarget> buffer_;
    std::size_t written_ = 0;
    bool closed_ = false;
  };

  /**
   * @brief Open a streaming session that writes to @p os.
   *
   * @param os          Destination stream.
   * @param threshold   Records buffered before an auto-flush
   *                    (default: DEFAULT_FLUSH_THRESHOLD).
   * @param write_header Emit the header immediately (default true).
   * @return A Session bound to @p os.
   */
  [[nodiscard]] Session
  open_session(std::ostream &os,
               std::size_t threshold = DEFAULT_FLUSH_THRESHOLD,
               bool write_header = true) const;

  /**
   * @brief Open a streaming session that writes to a new file at @p path.
   *
   * @param path      Destination file path (truncated).
   * @param threshold Records buffered before an auto-flush
   *                  (default: DEFAULT_FLUSH_THRESHOLD).
   * @return A Session that owns the file stream.
   * @throws std::runtime_error if the file cannot be opened.
   */
  [[nodiscard]] Session
  open_session_to_file(const std::string &path,
                       std::size_t threshold = DEFAULT_FLUSH_THRESHOLD) const;

private:
  std::unique_ptr<OffTargetFormatter> formatter_;
};

// =========================================================================
// ProfileWriter
// =========================================================================

/**
 * @brief Writes per-guide profiling statistics to the five legacy profile
 * files.
 *
 * @c ProfileWriter is stateless: it owns no data and performs no
 * accumulation.  It receives an already-built @c std::vector<GuideProfile>
 * (produced by @c ProfileAccumulator::build()) and writes the exact same
 * set of TSV-formatted @c .xls files that the legacy @c searchOnTST /
 * @c saveProfileGuide pipeline produced.
 *
 * ## Output files
 *
 * | Method                  | File suffix              | Description |
 * |-------------------------|--------------------------|-----------------------------------|
 * | write_profile()         | .profile.xls             | Per-position MM +
 * MM/OFFT summary | | write_extended_profile()| .extended_profile.xls    |
 * Per-(mm,nuc,pos) breakdown        | | write_profile_dna()     |
 * .profile_dna.xls         | DNA-bulge channel summary         | |
 * write_profile_rna()     | .profile_rna.xls         | RNA-bulge channel
 * summary         | | write_profile_complete()| .profile_complete.xls    |
 * All-channel combined summary      | | write_all_profiles()    | all five
 * files            | Convenience: writes all at once   |
 *
 * ## Column layouts
 *
 * **profile.xls** (one row per guide):
 * @code
 * GUIDE \t BP \t … \t BP \t \t ONT \t OFFT \t MM/OFFT \t \t 0MM \t 1MM \t …
 * @endcode
 * - BP columns (guide_len of them): @c pos_mm_count[i] — mismatch count at body
 * pos i.
 * - ONT: @c ont_count (0-MM, 0-bulge hits).
 * - OFFT: sum of @c offt_by_mm[1..max_mm].
 * - MM/OFFT: average mismatches per off-target (0 when OFFT == 0).
 * - nMM columns: @c offt_by_mm[0..max_mm].
 *
 * **extended_profile.xls** (one block per guide):
 * @code
 * >{guide+PAM_Ns}
 * \t BP \t … \t BP \t TARGETS
 * n MISMATCHES \t {ext_total[n][0..L-1]} \t {offt_complete[n]}
 * NUCLEOTIDE A \t {ext_nuc[n][0][0..L-1]}
 * NUCLEOTIDE C \t {ext_nuc[n][1][0..L-1]}
 * NUCLEOTIDE G \t {ext_nuc[n][2][0..L-1]}
 * NUCLEOTIDE T \t {ext_nuc[n][3][0..L-1]}
 * Bulge DNA    \t {ext_dna[n][0..L-1]}
 * Bulge RNA    \t {ext_rna[n][0..L-1]}
 * @endcode
 *
 * **profile_dna.xls / profile_rna.xls** (one row per guide):
 * @code
 * GUIDE \t BP \t … \t ONT \t OFFT \t MM/OFFT \t \t 0MM(1) \t 0MM(2) \t …
 * @endcode
 * - BP cells: @c pos_mm_in_dna[i](pos_bulge_dna[i]) format.
 * - nMM(b) columns: @c offt_dna[mm][b-1].
 *
 * **profile_complete.xls** (same shape as profile.xls but all-channel counts):
 * @code
 * GUIDE \t BP \t … \t ONT \t OFFT \t MM/OFFT \t \t 0MM \t 1MM \t …
 * @endcode
 *
 * ## Error handling
 * All methods throw @c std::runtime_error if a file cannot be opened or a
 * stream write fails, consistent with the rest of the output layer.
 *
 * ## Thread safety
 * @c ProfileWriter is stateless and reentrant.  Multiple threads may call
 * its methods concurrently as long as they write to different paths.
 */
class ProfileWriter {
public:
  // ---- stream API (testable without the filesystem) ------------------

  /**
   * @brief Write the main profile table (.profile.xls) to @p os.
   *
   * @param profiles  One @c GuideProfile per guide, in guide order.
   * @param os        Destination stream (must be good() on entry).
   * @throws std::runtime_error if the stream enters a fail state.
   */
  void write_profile(const std::vector<GuideProfile> &profiles,
                     std::ostream &os) const;

  /**
   * @brief Write the extended profile (.extended_profile.xls) to @p os.
   *
   * Produces one block per guide: a FASTA-style header line followed by
   * one sub-block per mismatch threshold (0..max_mm), each containing
   * seven rows (total, A, C, G, T, BulgeDNA, BulgeRNA).
   *
   * @param profiles  One @c GuideProfile per guide.
   * @param os        Destination stream.
   * @throws std::runtime_error if the stream enters a fail state.
   */
  void write_extended_profile(const std::vector<GuideProfile> &profiles,
                              std::ostream &os) const;

  /**
   * @brief Write the DNA-bulge channel profile (.profile_dna.xls) to @p os.
   *
   * BP cells are formatted as @c mm_count(bulge_count) matching the legacy
   * @c profiling_dna_mm[i](profiling_dna[i]) convention.
   *
   * @param profiles  One @c GuideProfile per guide.
   * @param os        Destination stream.
   * @throws std::runtime_error if the stream enters a fail state.
   */
  void write_profile_dna(const std::vector<GuideProfile> &profiles,
                         std::ostream &os) const;

  /**
   * @brief Write the RNA-bulge channel profile (.profile_rna.xls) to @p os.
   *
   * Same layout as write_profile_dna() but using RNA-channel counters.
   *
   * @param profiles  One @c GuideProfile per guide.
   * @param os        Destination stream.
   * @throws std::runtime_error if the stream enters a fail state.
   */
  void write_profile_rna(const std::vector<GuideProfile> &profiles,
                         std::ostream &os) const;

  /**
   * @brief Write the all-channel combined profile (.profile_complete.xls)
   *        to @p os.
   *
   * Same column layout as write_profile() but uses @c offt_complete_by_mm,
   * @c pos_mm_complete, and @c ont_count_complete, which aggregate across
   * MM-only, DNA-bulge, and RNA-bulge hits.
   *
   * @param profiles  One @c GuideProfile per guide.
   * @param os        Destination stream.
   * @throws std::runtime_error if the stream enters a fail state.
   */
  void write_profile_complete(const std::vector<GuideProfile> &profiles,
                              std::ostream &os) const;

  // ---- file API ------------------------------------------------------

  /**
   * @brief Write all five profile files in one call.
   *
   * Opens and writes:
   *   @c {stem}.profile.xls
   *   @c {stem}.extended_profile.xls
   *   @c {stem}.profile_dna.xls
   *   @c {stem}.profile_rna.xls
   *   @c {stem}.profile_complete.xls
   *
   * Each file is opened with @c std::ios::out | std::ios::trunc and closed
   * via RAII. The directory component of @p stem must already exist.
   *
   * @param profiles   One @c GuideProfile per guide, in guide order.
   * @param path_stem  Path prefix shared by all five output files.
   * @throws std::runtime_error if any file cannot be opened or written.
   */
  void write_all_profiles(const std::vector<GuideProfile> &profiles,
                          const std::string &path_stem) const;

private:
  // ---- internal helpers ----------------------------------------------

  /**
   * @brief Build the PAM-annotated guide string for header rows.
   *
   * Prepends or appends @c p.pam_len @c 'N' characters to @c p.guide
   * according to @c p.pam_at_start, matching the legacy guide-row format.
   *
   * @param p  The guide profile to annotate.
   * @return   Guide body with PAM placeholder Ns.
   * @complexity O(guide_len + pam_len).
   */
  [[nodiscard]] static std::string annotated_guide(const GuideProfile &p);

  /**
   * @brief Open a file or throw a descriptive runtime_error.
   *
   * @param path  File path to open (truncating).
   * @return      An open @c std::ofstream.
   * @throws std::runtime_error if the path cannot be opened.
   */
  [[nodiscard]] static std::ofstream open_or_throw(const std::string &path);

  /**
   * @brief Throw if @p os is in a fail or bad state.
   *
   * @param os     Stream to check.
   * @param where  Context label for the error message.
   * @throws std::runtime_error on stream failure.
   */
  static void check_stream(const std::ostream &os, const char *where);

  /**
   * @brief Write the common profile/complete header row to @p os.
   *
   * Emits: GUIDE \t BP×guide_len \t \t ONT \t OFFT \t MM/OFFT \t \t
   *        0MM \t 1MM \t … \t max_mm MM
   *
   * @param guide_len  Number of BP columns.
   * @param max_mm     Highest mismatch bucket label.
   * @param os         Destination stream.
   */
  static void write_mm_header(int guide_len, int max_mm, std::ostream &os);

  /**
   * @brief Write the DNA/RNA bulge profile header row to @p os.
   *
   * Emits: GUIDE \t BP×guide_len \t \t ONT \t OFFT \t MM/OFFT \t \t
   *        0MM(1) \t 0MM(2) \t … \t max_mm MM(max_bulge)
   *
   * @param guide_len  Number of BP columns.
   * @param max_mm     Highest mismatch bucket.
   * @param max_bulge  Highest bulge-size bucket.
   * @param os         Destination stream.
   */
  static void write_bulge_header(int guide_len, int max_mm, int max_bulge,
                                 std::ostream &os);
};

} // namespace crispritz