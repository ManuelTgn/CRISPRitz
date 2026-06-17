#pragma once

#include "offtarget.hpp"            // crispritz::OffTarget
#include "search_configuration.hpp"  // crispritz::OutputFormat
#include "tst_search.hpp"            // crispritz::SearchResult

#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace crispritz
{

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
    class OffTargetFormatter
    {
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
        [[nodiscard]] virtual std::string format_row(const OffTarget& ot) const = 0;

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
    class TsvFormatter final : public OffTargetFormatter
    {
      public:
        [[nodiscard]] std::string      header() const override;
        [[nodiscard]] std::string      format_row(const OffTarget& ot) const override;
        [[nodiscard]] std::string_view name() const noexcept override { return "tsv"; }
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
    class TargetsFormatter final : public OffTargetFormatter
    {
      public:
        [[nodiscard]] std::string      header() const override;
        [[nodiscard]] std::string      format_row(const OffTarget& ot) const override;
        [[nodiscard]] std::string_view name() const noexcept override { return "targets"; }
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
    [[nodiscard]] std::unique_ptr<OffTargetFormatter> make_formatter(OutputFormat fmt);

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
    class OutputWriter
    {
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
        std::size_t write(const std::vector<OffTarget>& records,
                          std::ostream&                 os,
                          bool                          write_header = true) const;

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
        std::size_t write(const SearchResult& result,
                          std::ostream&       os,
                          bool                write_header = true) const;

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
        std::size_t write_to_file(const std::vector<OffTarget>& records,
                                  const std::string&            path) const;

        /**
         * @brief Write a SearchResult to a file path.
         * @param result  The search result to serialize.
         * @param path    Destination file path.
         * @return        Total number of record rows written.
         * @throws std::runtime_error if the file cannot be opened or a write fails.
         */
        std::size_t write_to_file(const SearchResult& result,
                                  const std::string&  path) const;

        /** @return The canonical name of the active format. */
        [[nodiscard]] std::string_view format_name() const noexcept
        {
            return formatter_->name();
        }

      private:
        std::unique_ptr<OffTargetFormatter> formatter_;
    };

} // namespace crispritz